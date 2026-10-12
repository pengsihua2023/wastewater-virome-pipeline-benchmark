#!/usr/bin/env python3
"""Analyses added in response to peer review of the combined paper.

Three things the review asked for that the existing scripts did not produce.

(M1) All three cross-modality rank correlations, computed in one place.
     Only the viral-F1 correlation was previously computed (in
     17_combined_figures.py); the recall and precision values quoted in the
     paper were carried over from an analysis that predated the GOTTCHA2
     evidence-stream correction, and the recall value was wrong in sign.

(M2) The same correlations over defensible alternative analysis sets. Illumina
     nf-core/taxtriage and nf-core/taxprofiler are byte-identical, so including
     both as independent observations in a seven-point correlation is not
     defensible; this reports the de-duplicated set and the subset given reads
     on both modalities.

(M4) A cross-modality analysis restricted to the 53 viral members that the
     PacBio simulation did NOT raise to the depth floor. On this subset,
     designed abundance still maps to delivered reads on both arms, so it tests
     whether the failures and the rank instability survive removal of the floor
     artefact rather than being produced by it.

Writes 16_review_response_analyses.xlsx.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"

METRICS = ["recall_mean", "precision_mean", "f1_mean"]
LABEL = {"recall_mean": "viral recall", "precision_mean": "viral precision",
         "f1_mean": "viral F1"}

# Illumina taxtriage duplicates taxprofiler byte-for-byte (Methods, Table 2 fn a)
DUPLICATE = "TaxTriage"
# pipelines given reads on BOTH modalities (Table 2)
READS_BOTH = ["TaxProfiler", "TaxTriage", "CLARK", "GOTTCHA2"]


def correlations(piv, label, n_note=""):
    rows = []
    for m in METRICS:
        x, y = piv[m]["Illumina"], piv[m]["PacBio"]
        if x.isna().any() or y.isna().any() or len(x) < 3:
            rows.append(dict(analysis_set=label, metric=LABEL[m], n=len(x),
                             rho=np.nan, p=np.nan, note=n_note))
            continue
        r, p = stats.spearmanr(x, y)
        rows.append(dict(analysis_set=label, metric=LABEL[m], n=len(x),
                         rho=round(float(r), 3), p=round(float(p), 3),
                         note=n_note))
    return rows


def main():
    agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    piv = agg.pivot_table(index="pipeline", columns="platform", values=METRICS)

    # ------------------------------------------------- M1 / M2 correlations
    rows = []
    rows += correlations(piv, "all 7 pipelines (as originally reported)",
                         "includes the byte-identical Illumina taxtriage/taxprofiler pair")
    rows += correlations(piv.drop(DUPLICATE),
                         "6 pipelines, Illumina duplicate removed",
                         "nf-core/taxtriage dropped; identical to taxprofiler on Illumina")
    rows += correlations(piv.loc[READS_BOTH],
                         "4 pipelines given reads on both modalities",
                         "excludes contig-fed arms; still contains the duplicate pair")
    rows += correlations(piv.loc[[p for p in READS_BOTH if p != DUPLICATE]],
                         "3 pipelines, reads both sides and de-duplicated",
                         "smallest defensible set")
    corr = pd.DataFrame(rows)

    # --------------------------------------- M4 non-floored viral subset
    au = pd.read_excel(AN / "10_simulation_config_audit.xlsx", "per_genome")
    v = au[au.category.str.startswith("Virus (")].copy()
    v["t"] = pd.to_numeric(v["taxid"], errors="coerce")
    assert len(v) == 151, len(v)
    unfloored = set(v.loc[~v.floored, "t"].dropna().astype(int))
    floored = set(v.loc[v.floored, "t"].dropna().astype(int))

    det = pd.read_csv(INT / "viral_detection_detail.tsv", sep="\t")
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gv = gt[gt.category.str.startswith("Virus (")]
    k2t = dict(zip(gv.match_key, pd.to_numeric(gv.taxid, errors="coerce")))

    def as_taxid(k):
        try:
            return int(float(k))
        except (TypeError, ValueError):
            x = k2t.get(str(k))
            return int(x) if pd.notna(x) else None

    tp = det[det.status == "TP"].copy()
    tp["t"] = tp["key"].map(as_taxid)

    sub_rows = []
    for (plat, pipe, samp), g in tp.groupby(["platform", "pipeline", "sample"]):
        hit = set(g.t.dropna().astype(int))
        sub_rows.append(dict(
            platform=plat, pipeline=pipe, sample=samp,
            tp_unfloored=len(hit & unfloored), n_unfloored=len(unfloored),
            tp_floored=len(hit & floored), n_floored=len(floored)))
    # pipelines that detected nothing never appear in the TP table
    for plat in ["Illumina", "PacBio"]:
        for pipe in agg.pipeline.unique():
            for samp in ["sample_0", "sample_1", "sample_2"]:
                if not any(r["platform"] == plat and r["pipeline"] == pipe
                           and r["sample"] == samp for r in sub_rows):
                    sub_rows.append(dict(platform=plat, pipeline=pipe, sample=samp,
                                         tp_unfloored=0, n_unfloored=len(unfloored),
                                         tp_floored=0, n_floored=len(floored)))
    sub = pd.DataFrame(sub_rows)
    sub["recall_unfloored"] = sub.tp_unfloored / sub.n_unfloored
    sub["recall_floored"] = sub.tp_floored / sub.n_floored
    subagg = (sub.groupby(["platform", "pipeline"])
              [["recall_unfloored", "recall_floored"]].mean().reset_index())

    up = subagg.pivot_table(index="pipeline", columns="platform",
                            values="recall_unfloored")
    r_all, p_all = stats.spearmanr(up["Illumina"], up["PacBio"])
    up6 = up.drop(DUPLICATE)
    r_6, p_6 = stats.spearmanr(up6["Illumina"], up6["PacBio"])
    unf_corr = pd.DataFrame([
        dict(analysis_set="unfloored 53 viruses, all 7", metric="viral recall",
             n=len(up), rho=round(float(r_all), 3), p=round(float(p_all), 3)),
        dict(analysis_set="unfloored 53 viruses, duplicate removed",
             metric="viral recall", n=len(up6), rho=round(float(r_6), 3),
             p=round(float(p_6), 3))])

    # ---- second-round review: is the floor split confounded with abundance?
    # The unfloored members are the high-abundance stratum, so the split above
    # cannot separate floor status from abundance. This repeats the comparison
    # inside the one designed-abundance bin that contains both strata.
    band = v[(v.relative_abundance_num >= 1e-3) & (v.relative_abundance_num < 1e-2)]
    b_un = set(band.loc[~band.floored, "t"].dropna().astype(int))
    b_fl = set(band.loc[band.floored, "t"].dropna().astype(int))
    brows = []
    for (plat, pipe, samp), g in tp.groupby(["platform", "pipeline", "sample"]):
        hit = set(g.t.dropna().astype(int))
        brows.append(dict(platform=plat, pipeline=pipe, sample=samp,
                          recall_unfloored=len(hit & b_un) / len(b_un),
                          recall_floored=len(hit & b_fl) / len(b_fl)))
    for plat in ["Illumina", "PacBio"]:
        for pipe in agg.pipeline.unique():
            for samp in ["sample_0", "sample_1", "sample_2"]:
                if not any(r["platform"] == plat and r["pipeline"] == pipe
                           and r["sample"] == samp for r in brows):
                    brows.append(dict(platform=plat, pipeline=pipe, sample=samp,
                                      recall_unfloored=0.0, recall_floored=0.0))
    matched = (pd.DataFrame(brows).groupby(["platform", "pipeline"])
               [["recall_unfloored", "recall_floored"]].mean().reset_index())
    matched["n_unfloored"], matched["n_floored"] = len(b_un), len(b_fl)

    # abundance imbalance between the two strata, for disclosure in the text
    imbalance = pd.DataFrame([
        dict(stratum="unfloored", n=len(v[~v.floored]),
             median_abundance_pct=v.loc[~v.floored, "relative_abundance_num"].median(),
             n_ge_1e_2=int((v.loc[~v.floored, "relative_abundance_num"] >= 1e-2).sum()),
             n_lt_1e_3=int((v.loc[~v.floored, "relative_abundance_num"] < 1e-3).sum())),
        dict(stratum="floored", n=len(v[v.floored]),
             median_abundance_pct=v.loc[v.floored, "relative_abundance_num"].median(),
             n_ge_1e_2=int((v.loc[v.floored, "relative_abundance_num"] >= 1e-2).sum()),
             n_lt_1e_3=int((v.loc[v.floored, "relative_abundance_num"] < 1e-3).sum()))])
    u_stat, u_p = stats.mannwhitneyu(v.loc[~v.floored, "relative_abundance_num"],
                                     v.loc[v.floored, "relative_abundance_num"])
    imbalance["mannwhitney_p"] = f"{u_p:.2e}"

    # SD across the three samples, for Table 9
    sd = (sub.groupby(["platform", "pipeline"])
          [["recall_unfloored", "recall_floored"]].std().reset_index()
          .rename(columns={"recall_unfloored": "sd_unfloored",
                           "recall_floored": "sd_floored"}))
    subagg = subagg.merge(sd, on=["platform", "pipeline"])

    with pd.ExcelWriter(AN / "16_review_response_analyses.xlsx") as xw:
        matched.to_excel(xw, sheet_name="abundance_matched_floor", index=False)
        imbalance.to_excel(xw, sheet_name="floor_abundance_imbalance", index=False)
        corr.to_excel(xw, sheet_name="rank_correlations", index=False)
        subagg.to_excel(xw, sheet_name="recall_floored_vs_not", index=False)
        unf_corr.to_excel(xw, sheet_name="unfloored_correlation", index=False)
        sub.to_excel(xw, sheet_name="per_sample_floor_split", index=False)

    corr.to_csv(INT / "rank_correlations.tsv", sep="\t", index=False)
    subagg.to_csv(INT / "recall_floor_split.tsv", sep="\t", index=False)

    print("=== M1/M2 cross-modality rank correlations")
    print(corr.to_string(index=False))
    print(f"\n=== M4 depth floor: {len(unfloored)} of 151 viral members NOT floored, "
          f"{len(floored)} floored")
    print(subagg.round(3).to_string(index=False))
    print("\n=== M4 rank correlation on the unfloored subset only")
    print(unf_corr.to_string(index=False))
    print("\n=== R2: the two strata differ sharply in designed abundance")
    print(imbalance.to_string(index=False))
    print("\n=== R2/Q1: floor comparison inside the 1e-3 to 1e-2 percent band "
          f"({len(b_un)} unfloored vs {len(b_fl)} floored), PacBio")
    print(matched[matched.platform == "PacBio"].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
