#!/usr/bin/env python3
"""Score every pipeline against the 708-member CAMISIM ground truth.

Primary matching rule (strict, reproducible, species-level):
  * A ground-truth member is a TRUE POSITIVE for a given (platform, pipeline,
    sample) if the pipeline reports its NCBI taxid at a species-or-below node
    carrying more than `threshold` clade reads.
  * A reported species-or-below node whose taxid is NOT in the ground truth is
    a FALSE POSITIVE.
  * Ground-truth members with no such node are FALSE NEGATIVES.
  * True negatives are undefined (the negative space is each tool's entire
    reference database), therefore specificity, accuracy and MCC are NOT
    computed. See 08_statistical_analysis.md.

sourmash reports GTDB lineage strings and no NCBI taxids; it is scored by
normalised species name and flagged with match_mode = "name".

Each pipeline contributes exactly ONE primary evidence stream per platform
(`PRIMARY`); other streams are kept for the sensitivity analysis.

Outputs (interim/):
  sample_level_metrics.tsv   metrics at every read threshold
  detection_detail.tsv       per-item TP/FN/FP at the primary threshold
  sensitivity_streams.tsv    metrics for every alternative evidence stream
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
INT = ROOT / "analysis_708" / "interim"

SPECIES_RANKS = {"S", "S1", "S2", "S3", "S4", "S5", "species", "strain", "serotype",
                 "no rank", "subspecies", "biotype", "isolate", "genotype"}
SAMPLES = ["sample_0", "sample_1", "sample_2"]
THRESHOLDS = [0, 9, 99]          # node kept when reads > threshold  (>=1, >=10, >=100)
PRIMARY_THRESHOLD = 0

PRIMARY = {
    ("Illumina", "CLARK"): "clark_full",
    ("Illumina", "GOTTCHA2"): "gottcha2_tsv",
    ("Illumina", "KrakenMetaReads-nf"): "kmr_merged",
    ("Illumina", "mag-nf"): "mag_reads",
    ("Illumina", "TaxProfiler"): "taxprofiler_kraken2",
    ("Illumina", "TaxTriage"): "taxtriage_kraken2",
    ("Illumina", "sourmash"): "sourmash_combined",
    ("PacBio", "CLARK"): "clark_full",
    ("PacBio", "GOTTCHA2"): "gottcha2_tsv",
    ("PacBio", "KrakenMetaReads-nf"): "kmr_flye",
    # NOTE (2026-09, regenerated data): the three PacBio nf-core/mag READ-level
    # reports used to be byte-identical across samples (a data bug), which is
    # why "mag_reads_SUSPECT" was kept out of PRIMARY and "mag_contigs" used
    # instead. That bug is fixed in the regenerated data (MD5-verified distinct
    # per-sample reports), so the read-level stream is primary here too,
    # matching 14_viral_only_scoring.py's PRIMARY_FILES choice.
    ("PacBio", "mag-nf"): "mag_reads_SUSPECT",
    ("PacBio", "TaxProfiler"): "taxprofiler_kraken2",
    ("PacBio", "TaxTriage"): "taxtriage_kraken2",
    ("PacBio", "sourmash"): "sourmash_combined",
}

# streams retained only as sensitivity analyses
SENSITIVITY = [
    ("Illumina", "CLARK", "clark_g0.2"),
    ("Illumina", "KrakenMetaReads-nf", "kmr_spades"),
    ("Illumina", "KrakenMetaReads-nf", "kmr_megahit"),
    ("Illumina", "TaxProfiler", "taxprofiler_bracken"),
    ("Illumina", "TaxProfiler", "taxprofiler_abundance"),
    ("PacBio", "CLARK", "clark_g0.2"),
    ("PacBio", "TaxProfiler", "taxprofiler_abundance"),
    # kept as a sensitivity stream: the pre-fix contig-level evidence, for
    # comparison against the now-primary read-level stream above.
    ("PacBio", "mag-nf", "mag_contigs"),
]


def norm_key(n):
    s = re.sub(r"^[a-z]__", "", str(n).strip())   # GTDB rank prefix used by sourmash
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def score_one(sub, gt_taxids, gt_keys, uses_taxid):
    if uses_taxid:
        observed = set(sub["taxid_num"].dropna().astype(int))
        gt = gt_taxids
        reads_by = sub.groupby(sub["taxid_num"]).reads.max()
    else:
        observed = set(sub["name_key"]) - {""}
        gt = gt_keys
        reads_by = sub.groupby("name_key").reads.max()
    tp, fp, fn = observed & gt, observed - gt, gt - observed
    TP, FP, FN = len(tp), len(fp), len(fn)
    recall = TP / (TP + FN) if TP + FN else np.nan
    prec = TP / (TP + FP) if TP + FP else np.nan
    f1 = (2 * prec * recall / (prec + recall)) if (prec and recall) else 0.0
    return dict(n_reported_species_nodes=len(observed), TP=TP, FP=FP, FN=FN,
                recall=recall, precision=prec, f1=f1,
                fdr=(FP / (TP + FP) if TP + FP else np.nan),
                fnr=(FN / (TP + FN) if TP + FN else np.nan)), tp, fp, fn, reads_by


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gt["taxid_num"] = pd.to_numeric(gt["taxid"], errors="coerce").astype("Int64")
    assert len(gt) == 708, len(gt)
    gt_taxids = set(gt["taxid_num"].dropna().astype(int))
    gt_keys = set(gt["match_key"])
    tax2row = {int(r.taxid_num): r for r in gt.itertuples() if pd.notna(r.taxid_num)}
    key2row = {r.match_key: r for r in gt.itertuples()}

    det = pd.read_csv(INT / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    det["taxid_num"] = pd.to_numeric(det["taxid"], errors="coerce").astype("Int64")
    det["reads"] = pd.to_numeric(det["reads"], errors="coerce").fillna(0.0)
    det["name_key"] = det["name"].map(norm_key)
    sm = det["parser"].isin(["sourmash_prok", "sourmash_viral"])
    det["parser_primary"] = det["parser"].where(~sm, "sourmash_combined")
    det["is_species_node"] = det["rank_reported"].astype(str).isin(SPECIES_RANKS)

    records, detail, sens = [], [], []

    def run(plat, pipe, parser, sink, keep_detail):
        pool = det[(det.platform == plat) & (det.pipeline == pipe) &
                   (det.parser_primary == parser) & det.is_species_node]
        if pool.empty:
            print(f"WARNING: no rows for {plat}/{pipe}/{parser}")
            return
        uses_taxid = bool(pool["taxid_num"].notna().any())
        for sample in SAMPLES:
            base = pool[pool["sample"] == sample]
            for thr in THRESHOLDS:
                sub = base[base["reads"] > thr]
                m, tp, fp, fn, reads_by = score_one(sub, gt_taxids, gt_keys, uses_taxid)
                sink.append(dict(platform=plat, pipeline=pipe, sample=sample,
                                 evidence_stream=parser, match_mode=
                                 "taxid" if uses_taxid else "name",
                                 min_reads=thr + 1, n_ground_truth=708, **m))
                if keep_detail and thr == PRIMARY_THRESHOLD:
                    for st, items in (("TP", tp), ("FN", fn), ("FP", fp)):
                        for i in sorted(items):
                            row = (tax2row.get(i) if uses_taxid else key2row.get(i))
                            detail.append(dict(
                                platform=plat, pipeline=pipe, sample=sample, status=st,
                                key=i,
                                name=row.std_name if row is not None else "",
                                category=row.category if row is not None else "not in ground truth",
                                design_reads=row.design_reads_num if row is not None else np.nan,
                                rel_abundance_pct=row.relative_abundance_num if row is not None else np.nan,
                                observed_reads=float(reads_by.get(i, np.nan))))

    for (plat, pipe), parser in PRIMARY.items():
        run(plat, pipe, parser, records, True)
    for plat, pipe, parser in SENSITIVITY:
        run(plat, pipe, parser, sens, False)

    met = pd.DataFrame(records)
    met.to_csv(INT / "sample_level_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(detail).to_csv(INT / "detection_detail.tsv", sep="\t", index=False)
    pd.DataFrame(sens).to_csv(INT / "sensitivity_streams.tsv", sep="\t", index=False)

    show = met[met.min_reads == PRIMARY_THRESHOLD + 1]
    print(show.drop(columns=["n_ground_truth", "min_reads"]).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
