#!/usr/bin/env python3
"""Aggregate metrics, domain-stratified metrics, abundance analysis, error patterns.

Writes the 05/06/07 deliverable workbooks and interim tables used by the
figure script.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"
SUP = AN / "_superseded"          # whole-community scoring, superseded
SUP.mkdir(parents=True, exist_ok=True)

PIPELINES = ["GOTTCHA2", "mag-nf", "TaxProfiler", "TaxTriage",
             "KrakenMetaReads-nf", "CLARK", "sourmash"]


def boot_ci(vals, n=10000, seed=42):
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=float)
    if len(v) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    means = rng.choice(v, size=(n, len(v)), replace=True).mean(axis=1)
    return tuple(np.percentile(means, [2.5, 97.5]))


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    met = pd.read_csv(INT / "sample_level_metrics.tsv", sep="\t")
    det = pd.read_csv(INT / "detection_detail.tsv", sep="\t", low_memory=False)

    prim = met[met.min_reads == 1].copy()

    # ---------------------------------------------------- aggregate by platform
    agg = (prim.groupby(["platform", "pipeline", "evidence_stream", "match_mode"])
           .agg(TP_mean=("TP", "mean"), TP_sd=("TP", "std"),
                FP_mean=("FP", "mean"), FP_sd=("FP", "std"),
                FN_mean=("FN", "mean"), FN_sd=("FN", "std"),
                recall_mean=("recall", "mean"), recall_sd=("recall", "std"),
                precision_mean=("precision", "mean"), precision_sd=("precision", "std"),
                f1_mean=("f1", "mean"), f1_sd=("f1", "std"), n_samples=("f1", "size"))
           .reset_index())
    for m in ["recall", "precision", "f1"]:
        cis = []
        for _, r in agg.iterrows():
            vals = prim[(prim.platform == r.platform) & (prim.pipeline == r.pipeline)][m]
            cis.append(boot_ci(vals))
        agg[f"{m}_ci_lo"] = [c[0] for c in cis]
        agg[f"{m}_ci_hi"] = [c[1] for c in cis]

    # ------------------------------------------------ platform paired contrast
    pivot = prim.pivot_table(index=["pipeline", "sample"], columns="platform",
                             values=["recall", "precision", "f1", "TP", "FP"])
    plat_rows = []
    for pipe in PIPELINES:
        sub = prim[prim.pipeline == pipe]
        il = sub[sub.platform == "Illumina"].set_index("sample")
        pb = sub[sub.platform == "PacBio"].set_index("sample")
        common = sorted(set(il.index) & set(pb.index))
        if not common:
            continue
        for m in ["recall", "precision", "f1", "TP", "FP"]:
            d = pb.loc[common, m].values - il.loc[common, m].values
            lo, hi = boot_ci(d)
            plat_rows.append(dict(pipeline=pipe, metric=m,
                                  illumina_mean=il.loc[common, m].mean(),
                                  pacbio_mean=pb.loc[common, m].mean(),
                                  diff_pacbio_minus_illumina=d.mean(),
                                  diff_sd=d.std(ddof=1), diff_ci_lo=lo, diff_ci_hi=hi,
                                  n_pairs=len(common)))
    plat = pd.DataFrame(plat_rows)

    # ------------------------------------------------------- domain stratified
    cat2domain = {
        "Bacteria/Archaea": "Bacteria/Archaea",
        "Bacteriophage": "Bacteriophage",
        "Virus (automated)": "Virus (non-phage)",
        "Virus (manual, no-rank)": "Virus (non-phage)",
        "Virus (manual, low-read)": "Virus (non-phage)",
        "Virus (manual, not detected)": "Virus (non-phage)",
    }
    gt["domain"] = gt["category"].map(cat2domain)
    dom_sizes = gt["domain"].value_counts().to_dict()
    d2 = det[det.status.isin(["TP", "FN"])].copy()
    d2["domain"] = d2["category"].map(cat2domain)
    dom = (d2.groupby(["platform", "pipeline", "sample", "domain", "status"])
           .size().unstack("status").fillna(0).reset_index())
    dom["n_truth"] = dom["domain"].map(dom_sizes).astype(float)
    dom["recall"] = dom["TP"] / dom["n_truth"]
    dom_agg = (dom.groupby(["platform", "pipeline", "domain"])
               .agg(recall_mean=("recall", "mean"), recall_sd=("recall", "std"),
                    TP_mean=("TP", "mean"), n_truth=("n_truth", "first")).reset_index())

    # ------------------------------------------------- abundance-stratified recall
    # Design abundance is expressed on a single community-wide scale by the
    # "Relative abundance in simulation (%)" column; design_reads are block-
    # scaled and are NOT comparable between the bacterial and viral sections.
    bins = [0, 1e-4, 1e-3, 1e-2, 1e-1, np.inf]
    labels = ["<1e-4%", "1e-4-1e-3%", "1e-3-1e-2%", "1e-2-0.1%", ">=0.1%"]
    gt["abund_bin"] = pd.cut(gt["relative_abundance_num"], bins=bins, right=False, labels=labels)
    d2["abund_bin"] = pd.cut(d2["rel_abundance_pct"], bins=bins, right=False, labels=labels)
    bin_sizes = gt["abund_bin"].value_counts().to_dict()
    ab = (d2.groupby(["platform", "pipeline", "sample", "abund_bin", "status"],
                     observed=True).size().unstack("status").fillna(0).reset_index())
    ab["n_truth"] = ab["abund_bin"].astype(str).map({str(k): v for k, v in bin_sizes.items()})
    ab["recall"] = ab["TP"] / ab["n_truth"]
    ab_agg = (ab.groupby(["platform", "pipeline", "abund_bin"], observed=True)
              .agg(recall_mean=("recall", "mean"), recall_sd=("recall", "std"),
                   n_truth=("n_truth", "first")).reset_index())

    # ---------------------------------------------- observed vs designed reads
    tp = det[(det.status == "TP") & det.observed_reads.notna() &
             det.rel_abundance_pct.notna() & (det.rel_abundance_pct > 0)].copy()
    cors = []
    for (c_plat, c_pipe, c_samp), g in tp.groupby(["platform", "pipeline", "sample"]):
        g = g[g.observed_reads > 0]
        if len(g) < 10:
            continue
        x, y = np.log10(g.rel_abundance_pct), np.log10(g.observed_reads)
        sp = stats.spearmanr(g.rel_abundance_pct, g.observed_reads)
        pe = stats.pearsonr(x, y)
        cors.append(dict(platform=c_plat, pipeline=c_pipe, sample=c_samp, n_tp_with_reads=len(g),
                         spearman_rho=sp.statistic, spearman_p=sp.pvalue,
                         pearson_r_log10=pe.statistic, pearson_p_log10=pe.pvalue,
                         median_log10_ratio=float(np.median(y - x))))
    cor = pd.DataFrame(cors)

    # ---------------------------------------------------------- error patterns
    n_runs = det[["platform", "pipeline", "sample"]].drop_duplicates().shape[0]
    fn = det[det.status == "FN"]
    miss = (fn.groupby(["key", "name", "category", "rel_abundance_pct"]).size()
            .reset_index(name="n_runs_missed"))
    miss["frac_runs_missed"] = miss["n_runs_missed"] / n_runs
    miss = miss.sort_values(["n_runs_missed", "rel_abundance_pct"], ascending=[False, False])

    tpd = det[det.status == "TP"]
    det_by = (tpd.groupby(["key", "name", "category", "rel_abundance_pct"])["pipeline"]
              .nunique().reset_index(name="n_pipelines_detecting"))
    uniq = det_by[det_by.n_pipelines_detecting == 1].merge(
        tpd.groupby("key")["pipeline"].agg(lambda s: sorted(set(s))[0]).rename("only_pipeline"),
        on="key", how="left")

    # platform-exclusive ground-truth detections
    pres = (tpd.groupby(["key", "name", "category", "rel_abundance_pct", "platform"])
            .size().unstack("platform").fillna(0))
    pres["illumina_only"] = (pres.get("Illumina", 0) > 0) & (pres.get("PacBio", 0) == 0)
    pres["pacbio_only"] = (pres.get("PacBio", 0) > 0) & (pres.get("Illumina", 0) == 0)
    plat_excl = pres.reset_index()

    # recurrent false positives (taxid-based pipelines only)
    fp = det[(det.status == "FP")]
    fp_rec = (fp.groupby(["platform", "pipeline", "key"]).size()
              .reset_index(name="n_samples").query("n_samples == 3"))
    fp_top = (fp.sort_values("observed_reads", ascending=False)
              .groupby(["platform", "pipeline"]).head(20))

    # ------------------------------------------------------------------ write
    with pd.ExcelWriter(SUP / "05_sample_level_metrics.xlsx") as xw:
        prim.to_excel(xw, sheet_name="primary_min1_read", index=False)
        met.to_excel(xw, sheet_name="all_read_thresholds", index=False)
        pd.read_csv(INT / "sensitivity_streams.tsv", sep="\t").to_excel(
            xw, sheet_name="alternative_streams", index=False)
        dom.to_excel(xw, sheet_name="domain_stratified_by_sample", index=False)
        ab.to_excel(xw, sheet_name="abundance_stratified_by_sample", index=False)

    with pd.ExcelWriter(SUP / "06_aggregate_metrics.xlsx") as xw:
        agg.to_excel(xw, sheet_name="per_platform_pipeline", index=False)
        plat.to_excel(xw, sheet_name="platform_contrast_paired", index=False)
        dom_agg.to_excel(xw, sheet_name="domain_recall", index=False)
        ab_agg.to_excel(xw, sheet_name="abundance_recall", index=False)
        cor.to_excel(xw, sheet_name="abundance_correlation", index=False)

    with pd.ExcelWriter(SUP / "07_error_analysis.xlsx") as xw:
        miss.to_excel(xw, sheet_name="missed_by_n_runs", index=False)
        miss[miss.n_runs_missed == n_runs].to_excel(xw, sheet_name="missed_by_all_runs", index=False)
        uniq.to_excel(xw, sheet_name="detected_by_one_pipeline", index=False)
        plat_excl.to_excel(xw, sheet_name="platform_exclusive", index=False)
        fp_rec.to_excel(xw, sheet_name="recurrent_false_positives", index=False)
        fp_top.to_excel(xw, sheet_name="top_false_positives_by_reads", index=False)

    agg.to_csv(INT / "aggregate_metrics.tsv", sep="\t", index=False)
    plat.to_csv(INT / "platform_contrast.tsv", sep="\t", index=False)
    dom_agg.to_csv(INT / "domain_recall.tsv", sep="\t", index=False)
    ab_agg.to_csv(INT / "abundance_recall.tsv", sep="\t", index=False)
    cor.to_csv(INT / "abundance_correlation.tsv", sep="\t", index=False)
    miss.to_csv(INT / "missed_summary.tsv", sep="\t", index=False)

    print("n runs:", n_runs)
    print(agg[["platform", "pipeline", "recall_mean", "precision_mean", "f1_mean"]]
          .to_string(index=False))
    print("\nmissed by ALL runs:", int((miss.n_runs_missed == n_runs).sum()))
    print(dom_agg.to_string(index=False))
    print(ab_agg.head(40).to_string(index=False))
    print(cor.groupby(["platform", "pipeline"]).spearman_rho.mean().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
