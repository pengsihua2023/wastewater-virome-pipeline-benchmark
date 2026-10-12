#!/usr/bin/env python3
"""Viral-only stratified tables used by the figure scripts.

Produces, for the 151 non-phage viral ground truth:
  interim/viral_abundance_recall.tsv     recall per designed-abundance bin
  interim/viral_abundance_correlation.tsv Spearman rho of observed reads vs
                                          designed abundance over viral TPs
  interim/viral_detection_overlap.tsv     per-virus detection by modality

Bins contain only viral members, so the denominators differ from the
whole-community binning used earlier in the project.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
INT = ROOT / "analysis_708" / "interim"

BINS = [0, 1e-4, 1e-3, 1e-2, np.inf]
LABELS = ["<1e-4%", "1e-4-1e-3%", "1e-3-1e-2%", ">=1e-2%"]
SPECIES_RANKS = {"S", "S1", "S2", "S3", "S4", "S5", "species", "strain", "serotype",
                 "no rank", "subspecies", "biotype", "isolate", "genotype"}


def norm_key(n):
    s = re.sub(r"^[a-z]__", "", str(n).strip())
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    v = gt[gt.category.str.startswith("Virus (")].copy()
    assert len(v) == 151, len(v)
    v["t"] = pd.to_numeric(v["taxid"], errors="coerce")
    v["bin"] = pd.cut(v["relative_abundance_num"], BINS, right=False,
                      labels=LABELS).astype(str)
    bin_n = v["bin"].value_counts().to_dict()
    t2bin = dict(zip(v.t, v["bin"]))
    t2ab = dict(zip(v.t, v.relative_abundance_num))
    key2t = dict(zip(v.match_key, v.t))          # sourmash reports names, not taxids

    det = pd.read_csv(INT / "viral_detection_detail.tsv", sep="\t")

    def to_taxid(k):
        """detection_detail carries taxids for six pipelines and name keys for sourmash."""
        try:
            return int(float(k))
        except (TypeError, ValueError):
            return key2t.get(str(k))

    det["taxid"] = det["key"].map(to_taxid)
    det["bin"] = det["taxid"].map(t2bin)

    # ------------------------------------------------ abundance-binned recall
    tp = det[(det.status == "TP") & det["bin"].notna()]
    rows = tp.groupby(["platform", "pipeline", "sample", "bin"]).size().reset_index(name="TP")
    rows["n_truth"] = rows["bin"].map(bin_n)
    rows["recall"] = rows["TP"] / rows["n_truth"]
    # cells with zero detections do not appear in the groupby; add them back
    full = pd.MultiIndex.from_product(
        [sorted(det.platform.unique()), sorted(det.pipeline.unique()),
         sorted(det["sample"].unique()), LABELS],
        names=["platform", "pipeline", "sample", "bin"]).to_frame(index=False)
    rows = full.merge(rows, on=["platform", "pipeline", "sample", "bin"], how="left")
    rows["TP"] = rows["TP"].fillna(0)
    rows["n_truth"] = rows["bin"].map(bin_n)
    rows["recall"] = rows["TP"] / rows["n_truth"]
    rows.to_csv(INT / "viral_abundance_recall_per_sample.tsv", sep="\t", index=False)

    agg = (rows.groupby(["platform", "pipeline", "bin"])
           .agg(recall_mean=("recall", "mean"), recall_sd=("recall", "std"),
                n_truth=("n_truth", "first")).reset_index())
    agg.to_csv(INT / "viral_abundance_recall.tsv", sep="\t", index=False)

    # ------------------------------- observed reads vs designed viral abundance
    long = pd.read_csv(INT / "viral_detections_long.tsv", sep="\t", low_memory=False)
    long = long[long.is_viral & long["rank"].isin(SPECIES_RANKS) & (long.reads > 0)].copy()
    long["taxid_num"] = pd.to_numeric(long["taxid"], errors="coerce")
    long["t"] = long["taxid_num"]
    miss = long["t"].isna()
    long.loc[miss, "t"] = long.loc[miss, "name"].map(norm_key).map(key2t)
    long = long[long["t"].isin(set(v.t))]
    long["abund"] = long["t"].map(t2ab)

    cors = []
    for (plat, pipe, samp), g in long.groupby(["platform", "pipeline", "sample"]):
        g = g.groupby("t").agg(reads=("reads", "max"), abund=("abund", "first"))
        if len(g) < 10:
            cors.append(dict(platform=plat, pipeline=pipe, sample=samp,
                             n_tp=len(g), spearman_rho=np.nan, spearman_p=np.nan))
            continue
        sp = stats.spearmanr(g.abund, g.reads)
        cors.append(dict(platform=plat, pipeline=pipe, sample=samp, n_tp=len(g),
                         spearman_rho=sp.statistic, spearman_p=sp.pvalue))
    pd.DataFrame(cors).to_csv(INT / "viral_abundance_correlation.tsv", sep="\t", index=False)
    long[["platform", "pipeline", "sample", "t", "name", "reads", "abund"]].to_csv(
        INT / "viral_tp_reads.tsv", sep="\t", index=False)

    # ------------------------------------------------------- modality overlap
    byplat = tp.groupby("taxid").platform.agg(lambda s: frozenset(s))
    ov = pd.DataFrame({
        "category": ["both modalities", "PacBio only", "Illumina only", "neither"],
        "n_viruses": [
            sum(1 for x in byplat if len(x) == 2),
            sum(1 for x in byplat if x == frozenset({"PacBio"})),
            sum(1 for x in byplat if x == frozenset({"Illumina"})),
            151 - byplat.shape[0]]})
    ov.to_csv(INT / "viral_detection_overlap.tsv", sep="\t", index=False)

    print("bin sizes:", {k: bin_n[k] for k in LABELS})
    for plat in ["Illumina", "PacBio"]:
        print(f"\n--- {plat} viral recall by designed abundance")
        print(agg[agg.platform == plat].pivot_table(
            index="pipeline", columns="bin", values="recall_mean")
            .reindex(columns=LABELS).round(3).to_string())
    print("\noverlap:")
    print(ov.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
