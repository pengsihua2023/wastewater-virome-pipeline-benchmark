#!/usr/bin/env python3
"""Virus-focused re-scoring: restrict BOTH the ground truth and the false
positives to the viral space.

The earlier analysis scored all 708 community members together. Because 490 of
them (69%) are bacteria and archaea, that diluted every metric: a pipeline that
recovers all bacteria and no viruses still scored recall 0.69. It also counted
bacterial false positives against a pipeline's viral performance.

This script re-scores with the viral space only:
  * ground truth  = the 151 non-phage viruses (primary) or all 218 viral
    members including the 67 ssRNA phage (secondary)
  * false positive = a reported species-or-below node that lies inside the
    Viruses branch and is not in the viral ground truth. When the ground truth
    is the 151 non-phage viruses, the 67 phage members of the community are
    OUT OF SCOPE and are removed from the reported set rather than counted as
    false positives.
  * bacterial/archaeal detections and bacterial/archaeal false positives are
    excluded from every metric

Determining whether a reported node is viral, per tool:
  Kraken2 reports (taxprofiler, taxtriage, mag-nf, KrakenMetaReads-nf)
      indentation-tracked descent from the "Viruses" (taxid 10239) D node
  CLARK        Lineage column begins with "Viruses;"
  sourmash     lineage begins with "Viruses" (viral DB) / "d__Viruses"
  GOTTCHA2     no lineage in the scored TSV, so viral membership is resolved
               from the sibling .lineage.tsv file

Outputs:
  interim/viral_detections_long.tsv    every reported node with is_viral
  interim/viral_sample_metrics.tsv     virus-only metrics per run
  11_viral_only_metrics.xlsx           deliverable workbook
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

SAMPLES = ["sample_0", "sample_1", "sample_2"]
SPECIES_RANKS = {"S", "S1", "S2", "S3", "S4", "S5", "species", "strain", "serotype",
                 "no rank", "subspecies", "biotype", "isolate", "genotype"}
THRESHOLDS = [0, 9, 99]

# one primary evidence file per pipeline per platform
PRIMARY_FILES = {
    ("Illumina", "TaxProfiler"): "Illumina/TaxProfiler/{s}/sample_708_{i}_null_Viral_ref.kraken2.kraken2.report.txt",
    ("Illumina", "TaxTriage"): "Illumina/TaxTriage/{s}/Illumina_708_{i}.kraken2.report.txt",
    ("Illumina", "mag-nf"): "Illumina/mag-nf/{s}/Sample_708_illumina_{i}.kraken2_report.txt",
    ("PacBio", "TaxProfiler"): "Pacbio/TaxProfiler/{s}/sample_708_{i}_PacBio_null_Viral_ref.kraken2.kraken2.report.txt",
    ("PacBio", "TaxTriage"): "Pacbio/TaxTriage/{s}/PacBio_708_{i}.kraken2.report.txt",
    # NOTE (2026-09, regenerated data): the three PacBio nf-core/mag READ-level
    # reports used to be byte-identical across samples (a data bug), which is
    # why earlier analyses substituted the contig-level report instead. That
    # bug is fixed in the regenerated data -- the three per-sample read-level
    # reports are now genuinely distinct (MD5-verified) -- so this now uses the
    # read-level report, matching the read-level evidence used for every other
    # PacBio pipeline except the contig-native ones (Table 2).
    ("PacBio", "mag-nf"): "Pacbio/mag-nf/{s}/sample_708_pacbio_{i}_reads_kraken2_report.txt",
}
# kept for reference / sensitivity analysis: the old contig-level substitute
PACBIO_MAGNF_CONTIGS = "Pacbio/mag-nf/{s}/sample_708_pacbio_{i}_contigs_kraken2_report.txt"


def kraken_viral_nodes(path):
    """Walk a Kraken2 report, returning rows with an is_viral flag.

    The report is depth-ordered and indentation encodes depth, so everything
    more deeply indented than the Viruses (10239) node belongs to that clade
    until the indentation returns to its level or shallower.
    """
    rows = []
    viral_indent = None
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        _pct, clade, _direct, rank, taxid, name = parts[:6]
        indent = len(name) - len(name.lstrip())
        nm = name.strip()
        if taxid.strip() == "10239":
            viral_indent = indent
            is_viral = True
        elif viral_indent is not None and indent > viral_indent:
            is_viral = True
        else:
            if viral_indent is not None and indent <= viral_indent:
                viral_indent = None
            is_viral = False
        try:
            reads = float(clade)
        except ValueError:
            continue
        rows.append(dict(taxid=taxid.strip(), name=nm, rank=rank.strip(),
                         reads=reads, is_viral=is_viral))
    return rows


def gottcha2_rows(tsv, lineage):
    """GOTTCHA2: viral membership resolved from the sibling .lineage.tsv."""
    viral_taxids = set()
    if lineage.exists():
        for line in open(lineage):
            f = line.rstrip("\n").split("\t")
            # fields repeat as rank, taxid, name; superkingdom is first
            if len(f) > 3 and f[1] == "superkingdom" and f[3] == "Viruses":
                for k in range(2, len(f) - 1, 3):
                    viral_taxids.add(f[k].strip())
    rows = []
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    for _, r in df.iterrows():
        t = str(r["TAXID"]).strip()
        rows.append(dict(taxid=t, name=str(r["NAME"]).strip(),
                         rank=str(r["LEVEL"]).strip(),
                         reads=float(r["READ_COUNT"]),
                         is_viral=(t in viral_taxids or t == "10239")))
    return rows


def clark_rows(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if not r.get("Name"):
                continue
            lin = (r.get("Lineage") or "")
            try:
                reads = float(r.get("Count") or 0)
            except ValueError:
                reads = 0.0
            rows.append(dict(taxid=(r.get("TaxID") or "").strip(),
                             name=r["Name"].strip(), rank="species", reads=reads,
                             is_viral=lin.startswith("Viruses")))
    return rows


def sourmash_rows(paths):
    rows = []
    for p in paths:
        df = pd.read_csv(p, dtype=str)
        sub = df[df["rank"] == "species"]
        for _, r in sub.iterrows():
            lin = str(r["lineage"])
            parts = [x.strip() for x in lin.split(";")]
            genus = parts[-2] if len(parts) >= 2 else ""
            genus = re.sub(r"^[a-z]__", "", genus)
            rows.append(dict(taxid="", name=parts[-1],
                             rank="species",
                             reads=float(pd.to_numeric(r["f_weighted_at_rank"],
                                                       errors="coerce") or 0) * 1e6,
                             is_viral=lin.startswith("Viruses") or
                                      lin.startswith("d__Viruses"),
                             genus_from_lineage=genus))
    return rows


def norm_key(n):
    s = re.sub(r"^[a-z]__", "", str(n).strip())
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def collect(platform, pipeline, sample, i):
    """Return the primary evidence rows for one run."""
    key = (platform, pipeline)
    if key in PRIMARY_FILES:
        return kraken_viral_nodes(ROOT / PRIMARY_FILES[key].format(s=sample, i=i))
    base = ROOT / ("Illumina" if platform == "Illumina" else "Pacbio")
    if pipeline == "KrakenMetaReads-nf":
        if platform == "Illumina":
            # merged two-assembler CSV carries no lineage; use the SPAdes report,
            # which is the deeper of the two, plus the MEGAHIT report
            rows = []
            for p in sorted((base / pipeline / sample).glob("*_report.txt")):
                if "merged" in p.name:
                    continue
                rows += kraken_viral_nodes(p)
            return rows
        return kraken_viral_nodes(
            next((base / pipeline / sample).glob("*_flye_report.txt")))
    if pipeline == "GOTTCHA2":
        # GOTTCHA2 writes a filtered ".tsv" and an unfiltered ".full.tsv".
        # The other six pipelines are scored on unfiltered output, so the
        # unfiltered file is used here for a like-for-like comparison; the
        # filtered file is retained as a sensitivity stream.
        d = base / pipeline / sample
        full = sorted(d.glob("*.full.tsv"))[0]
        lineage = full.parent / full.name.replace(".full.tsv", ".lineage.tsv")
        return gottcha2_rows(full, lineage)
    if pipeline == "CLARK":
        d = base / pipeline / sample
        f = [p for p in sorted(d.glob("*abundance*.csv")) if ".g0.2." not in p.name][0]
        return clark_rows(f)
    if pipeline == "sourmash":
        d = base / pipeline / sample
        return sourmash_rows(sorted(d.glob("*.summarized.csv")))
    raise KeyError(pipeline)


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gt["taxid_num"] = pd.to_numeric(gt["taxid"], errors="coerce").astype("Int64")
    nonphage = gt[gt.category.str.startswith("Virus (")]
    phage = gt[gt.category == "Bacteriophage"]
    assert len(nonphage) == 151 and len(phage) == 67

    # When scoring the 151 non-phage viruses, the 67 phage members are part of
    # the simulated community but out of scope: they are neither true nor false
    # positives, so they are EXCLUDED from the reported nodes rather than
    # counted against the pipeline.
    phage_tax = set(phage.taxid_num.dropna().astype(int))
    phage_key = set(phage.match_key)
    allv = pd.concat([nonphage, phage])
    GT = {
        "non-phage virus (n=151)": (set(nonphage.taxid_num.dropna().astype(int)),
                                    set(nonphage.match_key),
                                    phage_tax, phage_key),
        "all virus incl. phage (n=218)": (
            set(allv.taxid_num.dropna().astype(int)),
            set(allv.match_key), set(), set()),
    }

    all_rows, records, detail = [], [], []
    PIPES = ["TaxProfiler", "TaxTriage", "mag-nf", "KrakenMetaReads-nf",
             "GOTTCHA2", "CLARK", "sourmash"]
    for platform in ["Illumina", "PacBio"]:
        for pipeline in PIPES:
            for i, sample in enumerate(SAMPLES):
                rows = collect(platform, pipeline, sample, i)
                df = pd.DataFrame(rows)
                df["platform"], df["pipeline"], df["sample"] = platform, pipeline, sample
                all_rows.append(df)

                uses_taxid = df["taxid"].astype(str).str.strip().ne("").any()
                df["taxid_num"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")
                df["name_key"] = df["name"].map(norm_key)
                # a ground-truth member is viral by definition, so recall is
                # evaluated over every species-or-below node; the is_viral flag
                # only bounds which unexpected nodes may count as false positives
                sp_all = df[df["rank"].isin(SPECIES_RANKS)]
                sp = sp_all[sp_all["is_viral"]]

                for gt_name, (tax_set, key_set, ex_tax, ex_key) in GT.items():
                    for thr in THRESHOLDS:
                        s_v = sp[sp["reads"] > thr]            # viral-branch nodes
                        s_a = sp_all[sp_all["reads"] > thr]    # all species nodes
                        if uses_taxid:
                            truth = tax_set
                            hit = set(s_a["taxid_num"].dropna().astype(int))
                            fp_pool = set(s_v["taxid_num"].dropna().astype(int))
                        else:
                            truth = key_set
                            hit = set(s_a["name_key"]) - {""}
                            fp_pool = set(s_v["name_key"]) - {""}
                        tp = hit & truth
                        fp = (fp_pool - truth) - (ex_tax if uses_taxid else ex_key)
                        fn = truth - tp
                        obs = tp | fp
                        TP, FP, FN = len(tp), len(fp), len(fn)
                        rec = TP / (TP + FN) if TP + FN else np.nan
                        pre = TP / (TP + FP) if TP + FP else np.nan
                        f1 = (2 * pre * rec / (pre + rec)) if (pre and rec) else 0.0
                        records.append(dict(
                            platform=platform, pipeline=pipeline, sample=sample,
                            ground_truth=gt_name, min_reads=thr + 1,
                            match_mode="taxid" if uses_taxid else "name",
                            n_truth=len(truth), n_viral_nodes_reported=len(obs),
                            TP=TP, FP=FP, FN=FN, recall=rec, precision=pre, f1=f1))
                        if thr == 0 and gt_name.startswith("non-phage"):
                            for st, items in (("TP", tp), ("FN", fn), ("FP", fp)):
                                for it in sorted(items):
                                    detail.append(dict(
                                        platform=platform, pipeline=pipeline,
                                        sample=sample, status=st, key=it))

    # ---- sensitivity: GOTTCHA2 scored on its FILTERED ".tsv" instead of the
    # unfiltered ".full.tsv" used above, so the effect of that choice is on record
    gsens = []
    for platform in ["Illumina", "PacBio"]:
        base = ROOT / ("Illumina" if platform == "Illumina" else "Pacbio")
        for sample in SAMPLES:
            d = base / "GOTTCHA2" / sample
            filt = [q for q in sorted(d.glob("*.tsv"))
                    if not q.name.endswith((".full.tsv", ".lineage.tsv"))][0]
            lin = filt.parent / filt.name.replace(".tsv", ".lineage.tsv")
            df = pd.DataFrame(gottcha2_rows(filt, lin))
            df["taxid_num"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")
            sp_all = df[df["rank"].isin(SPECIES_RANKS)]
            sp_v = sp_all[sp_all["is_viral"]]
            tax_set, _, ex_tax, _ = GT["non-phage virus (n=151)"]
            for thr in THRESHOLDS:
                hit = set(sp_all[sp_all.reads > thr]["taxid_num"].dropna().astype(int))
                fpp = set(sp_v[sp_v.reads > thr]["taxid_num"].dropna().astype(int))
                tp = hit & tax_set
                fp = (fpp - tax_set) - ex_tax
                TP, FP, FN = len(tp), len(fp), len(tax_set) - len(tp)
                rec = TP / (TP + FN) if TP + FN else np.nan
                pre = TP / (TP + FP) if TP + FP else np.nan
                gsens.append(dict(platform=platform, pipeline="GOTTCHA2 (filtered .tsv)",
                                  sample=sample, min_reads=thr + 1, TP=TP, FP=FP, FN=FN,
                                  recall=rec, precision=pre,
                                  f1=(2 * pre * rec / (pre + rec)) if (pre and rec) else 0.0))
    gsens = pd.DataFrame(gsens)
    gsens.to_csv(INT / "gottcha2_filtered_sensitivity.tsv", sep="\t", index=False)

    long = pd.concat(all_rows, ignore_index=True)
    long.to_csv(INT / "viral_detections_long.tsv", sep="\t", index=False)
    met = pd.DataFrame(records)
    met.to_csv(INT / "viral_sample_metrics.tsv", sep="\t", index=False)
    dt = pd.DataFrame(detail)
    dt.to_csv(INT / "viral_detection_detail.tsv", sep="\t", index=False)

    # ------------------------------------------------------------- aggregate
    prim = met[(met.min_reads == 1) &
               (met.ground_truth == "non-phage virus (n=151)")]
    agg = (prim.groupby(["platform", "pipeline"])
           .agg(TP_mean=("TP", "mean"), TP_sd=("TP", "std"),
                FP_mean=("FP", "mean"), FP_sd=("FP", "std"),
                FN_mean=("FN", "mean"),
                recall_mean=("recall", "mean"), recall_sd=("recall", "std"),
                precision_mean=("precision", "mean"), precision_sd=("precision", "std"),
                f1_mean=("f1", "mean"), f1_sd=("f1", "std"))
           .reset_index())
    agg.to_csv(INT / "viral_aggregate_metrics.tsv", sep="\t", index=False)

    # ------------------------------------------------- viral error analysis
    # named/annotated views of the per-run TP/FN/FP detail, so the viral
    # workbook is self-contained and does not depend on the whole-community
    # error analysis (which uses a different, superseded scoring scope)
    v = nonphage.copy()
    v["t"] = pd.to_numeric(v["taxid"], errors="coerce")
    k2t = dict(zip(v.match_key, v.t))
    t2name = dict(zip(v.t, v.std_name))
    t2ab = dict(zip(v.t, v.relative_abundance_num))
    t2cat = dict(zip(v.t, v.category))

    dd = dt.copy()

    def as_taxid(k):
        try:
            return int(float(k))
        except (TypeError, ValueError):
            return k2t.get(str(k))

    dd["taxid"] = dd["key"].map(as_taxid)
    dd["virus"] = dd["taxid"].map(t2name)
    dd["designed_abundance_pct"] = dd["taxid"].map(t2ab)
    dd["category"] = dd["taxid"].map(t2cat)

    n_runs_per_plat = 21          # 7 pipelines x 3 samples
    fn = dd[dd.status == "FN"]
    missed = (fn.groupby(["platform", "taxid", "virus", "designed_abundance_pct"])
              .size().reset_index(name="n_runs_missed"))
    missed["frac_runs_missed"] = missed["n_runs_missed"] / n_runs_per_plat
    missed = missed.sort_values(["platform", "n_runs_missed", "designed_abundance_pct"],
                                ascending=[True, False, False])

    tpd = dd[dd.status == "TP"]
    by_pipe = (tpd.groupby(["platform", "taxid", "virus", "designed_abundance_pct"])
               ["pipeline"].agg(lambda s: "; ".join(sorted(set(s)))).reset_index()
               .rename(columns={"pipeline": "detected_by"}))

    sentinel = by_pipe[by_pipe["taxid"].map(t2cat).astype(str)
                       .str.startswith("Virus (manual")]

    modality = (tpd.groupby("taxid").platform.agg(lambda s: frozenset(s))
                .reset_index(name="platforms"))
    modality["virus"] = modality["taxid"].map(t2name)
    modality["designed_abundance_pct"] = modality["taxid"].map(t2ab)
    modality["detected_on"] = modality["platforms"].map(
        lambda x: "both" if len(x) == 2 else next(iter(x)))
    never = sorted(set(v.t.dropna()) - set(tpd["taxid"].dropna()))
    never_df = pd.DataFrame({"taxid": never})
    never_df["virus"] = never_df["taxid"].map(t2name)
    never_df["designed_abundance_pct"] = never_df["taxid"].map(t2ab)

    with pd.ExcelWriter(AN / "11_viral_only_metrics.xlsx") as xw:
        agg.to_excel(xw, sheet_name="viral_aggregate_151", index=False)
        prim.to_excel(xw, sheet_name="per_sample_151_min1", index=False)
        met.to_excel(xw, sheet_name="all_thresholds_both_GT", index=False)
        (long.groupby(["platform", "pipeline", "sample"])
         .agg(nodes_total=("name", "size"),
              nodes_viral=("is_viral", "sum")).reset_index()
         .to_excel(xw, sheet_name="viral_node_counts", index=False))
        missed.to_excel(xw, sheet_name="viruses_missed_per_run", index=False)
        by_pipe.to_excel(xw, sheet_name="viruses_detected_by", index=False)
        sentinel.to_excel(xw, sheet_name="sentinel_virus_detection", index=False)
        modality[["taxid", "virus", "designed_abundance_pct", "detected_on"]].to_excel(
            xw, sheet_name="modality_overlap", index=False)
        never_df.to_excel(xw, sheet_name="never_detected", index=False)
        gsens.to_excel(xw, sheet_name="gottcha2_filtered_sensitivity", index=False)

    print("=== VIRUS-ONLY metrics, ground truth = 151 non-phage viruses, >=1 read")
    for plat in ["Illumina", "PacBio"]:
        a = agg[agg.platform == plat].sort_values("f1_mean", ascending=False)
        print(f"\n--- {plat}")
        for _, r in a.iterrows():
            print(f"  {r.pipeline:20s} TP={r.TP_mean:6.1f}  FP={r.FP_mean:8.1f}  "
                  f"FN={r.FN_mean:6.1f}  R={r.recall_mean:.4f}  "
                  f"P={r.precision_mean:.4f}  F1={r.f1_mean:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
