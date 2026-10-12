#!/usr/bin/env python3
"""Two audits that change how pipeline failures should be attributed.

(1) CLARK's confidence signal is gamma, not read count.
    CLARK writes an unfiltered abundance table and a gamma>=0.2 filtered one.
    The unfiltered table was used as the primary evidence for a like-for-like
    comparison with the other pipelines, but CLARK's own run scripts set
    MODE=0 specifically to obtain the gamma column and filter on it. This
    script scores both, and shows that a read-count threshold cannot substitute
    for the gamma filter because CLARK's true and false positives have
    overlapping read-count distributions.

(2) What each pipeline was actually given as input.
    sourmash's `summarized.csv` records the query file it sketched. On Illumina
    that is a reads FASTQ; on PacBio it is `metaflye/assembly.fasta`. Attributing
    sourmash's PacBio result to "long reads" is therefore wrong: it was given
    assembled contigs, like nf-core/mag and KrakenMetaReads-nf. Grouping the
    PacBio viral recall by input type rather than by pipeline family makes the
    pattern unambiguous.

Writes 14_clark_gamma_and_input_types.xlsx.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

SAMPLES = ["sample_0", "sample_1", "sample_2"]
THRESHOLDS = [0, 9, 99]


def clark_files(platform, sample):
    base = ROOT / ("Illumina" if platform == "Illumina" else "Pacbio") / "CLARK" / sample
    unfilt = [p for p in sorted(base.glob("*abundance*.csv")) if ".g0.2." not in p.name][0]
    filt = sorted(base.glob("*.g0.2.csv"))[0]
    return unfilt, filt


def score_clark(path, vt, ph):
    d = pd.read_csv(path)
    d["t"] = pd.to_numeric(d["TaxID"], errors="coerce")
    d["viral"] = d["Lineage"].astype(str).str.startswith("Viruses")
    out = {}
    for thr in THRESHOLDS:
        s = d[d.viral & (d["Count"] > thr)]
        obs = set(s.t.dropna().astype(int)) - ph
        tp, fp = obs & vt, obs - vt
        TP, FP = len(tp), len(fp)
        rec = TP / len(vt)
        pre = TP / (TP + FP) if TP + FP else np.nan
        out[thr + 1] = dict(TP=TP, FP=FP, recall=rec, precision=pre,
                            f1=(2 * pre * rec / (pre + rec)) if (pre and rec) else 0.0)
    return out, d


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    v = gt[gt.category.str.startswith("Virus (")]
    vt = set(pd.to_numeric(v["taxid"], errors="coerce").dropna().astype(int))
    ph = set(pd.to_numeric(gt[gt.category == "Bacteriophage"]["taxid"],
                           errors="coerce").dropna().astype(int))

    # ---------------------------------------------- (1) CLARK gamma filter
    rows, readdist = [], []
    for platform in ["Illumina", "PacBio"]:
        for sample in SAMPLES:
            unf, flt = clark_files(platform, sample)
            for label, path in [("unfiltered", unf), ("gamma>=0.2", flt)]:
                res, d = score_clark(path, vt, ph)
                for thr, m in res.items():
                    rows.append(dict(platform=platform, sample=sample,
                                     stream=label, min_reads=thr, **m))
                if label == "unfiltered":
                    vir = d[d.viral & (d["Count"] > 0)]
                    tpc = vir[vir.t.isin(vt)]["Count"]
                    fpc = vir[~vir.t.isin(vt) & ~vir.t.isin(ph)]["Count"]
                    readdist.append(dict(platform=platform, sample=sample,
                                         tp_n=len(tpc), tp_median=tpc.median(),
                                         fp_n=len(fpc), fp_median=fpc.median(),
                                         fp_frac_single_read=float((fpc == 1).mean())))
    clark = pd.DataFrame(rows)
    clark_agg = (clark.groupby(["platform", "stream", "min_reads"])
                 [["TP", "FP", "recall", "precision", "f1"]].mean().reset_index())
    rd = pd.DataFrame(readdist)

    # ---------------------------------------- (2) input type per pipeline
    inputs = []
    for platform, folder in [("Illumina", "Illumina"), ("PacBio", "Pacbio")]:
        for sample in SAMPLES:
            f = sorted((ROOT / folder / "sourmash" / sample).glob("*.viral.summarized.csv"))
            if f:
                q = pd.read_csv(f[0], nrows=1)["query_filename"].iloc[0]
                inputs.append(dict(platform=platform, pipeline="sourmash", sample=sample,
                                   query_file=Path(str(q)).name,
                                   input_type="contigs" if "assembly" in str(q) else "reads",
                                   source="query_filename field of summarized.csv"))
    # the remaining pipelines: input type read from their run scripts / evidence choice
    declared = {
        ("Illumina", "TaxProfiler"): ("reads", "read-level Kraken2 report"),
        ("Illumina", "TaxTriage"): ("reads", "read-level Kraken2 report"),
        ("Illumina", "mag-nf"): ("reads", "read-level Kraken2 report"),
        ("Illumina", "KrakenMetaReads-nf"): ("contigs", "metaSPAdes + MEGAHIT contigs"),
        ("Illumina", "GOTTCHA2"): ("reads", "-xm sr on FASTQ"),
        ("Illumina", "CLARK"): ("reads", "classify_metagenome.sh on FASTQ"),
        ("PacBio", "TaxProfiler"): ("reads", "read-level Kraken2 report"),
        ("PacBio", "TaxTriage"): ("reads", "read-level Kraken2 report"),
        ("PacBio", "mag-nf"): ("contigs", "contig-level Kraken2 report"),
        ("PacBio", "KrakenMetaReads-nf"): ("contigs", "metaFlye contigs"),
        ("PacBio", "GOTTCHA2"): ("reads", "-xm map-hifi on FASTQ"),
        ("PacBio", "CLARK"): ("reads", "classify_metagenome.sh on FASTQ"),
    }
    for (plat, pipe), (typ, src) in declared.items():
        inputs.append(dict(platform=plat, pipeline=pipe, sample="all",
                           query_file="", input_type=typ, source=src))
    inp = pd.DataFrame(inputs)

    itype = (inp[inp["sample"].isin(["all", "sample_0"])]
             .drop_duplicates(["platform", "pipeline"])[["platform", "pipeline", "input_type"]])
    agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    byinput = agg.merge(itype, on=["platform", "pipeline"])[
        ["platform", "input_type", "pipeline", "recall_mean", "precision_mean", "f1_mean"]]
    byinput = byinput.sort_values(["platform", "input_type", "recall_mean"],
                                  ascending=[True, True, False])

    with pd.ExcelWriter(AN / "14_clark_gamma_and_input_types.xlsx") as xw:
        clark_agg.to_excel(xw, sheet_name="clark_gamma_vs_unfiltered", index=False)
        clark.to_excel(xw, sheet_name="clark_per_sample", index=False)
        rd.to_excel(xw, sheet_name="clark_read_count_overlap", index=False)
        inp.to_excel(xw, sheet_name="pipeline_input_types", index=False)
        byinput.to_excel(xw, sheet_name="viral_recall_by_input_type", index=False)

    clark_agg.to_csv(INT / "clark_gamma_comparison.tsv", sep="\t", index=False)
    byinput.to_csv(INT / "viral_recall_by_input_type.tsv", sep="\t", index=False)

    print("=== CLARK: unfiltered vs gamma>=0.2 (mean of 3 samples, >=1 read)")
    print(clark_agg[clark_agg.min_reads == 1].round(4).to_string(index=False))
    print("\n=== CLARK: read counts cannot separate TP from FP (unfiltered)")
    print(rd.groupby("platform")[["tp_median", "fp_median", "fp_frac_single_read"]]
          .mean().round(3).to_string())
    print("\n=== sourmash input files")
    print(inp[inp.pipeline == "sourmash"][["platform", "sample", "query_file", "input_type"]]
          .to_string(index=False))
    print("\n=== PacBio viral recall grouped by input type")
    print(byinput[byinput.platform == "PacBio"].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
