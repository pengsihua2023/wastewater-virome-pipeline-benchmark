#!/usr/bin/env python3
"""Additional analyses requested in the pre-submission review (round 1), using existing data only.

1. bootstrap_ci.tsv             95% bootstrap CIs over the 151 viruses (recall, F1 with FP fixed at its observed mean) at
                                >=1 and >=10 reads, and paired modality differences (PacBio - Illumina) for each pipeline.
2. floor_excluded_metrics.tsv   metrics recomputed on the 98 viruses NOT affected by the one-read minimum.
3. false_positive_profile.tsv   read support and taxonomic relatedness of strict false positives (>=1 read).
4. kraken2_misses_illumina.tsv  viruses missed by nf-core/taxprofiler on Illumina, whether their taxid exists in the
                                Kraken2 database report, and which other pipelines detect them.

Inputs: interim/coverage/per_virus_detection_by_threshold.tsv (samples detecting each virus, 0-3),
interim/viral_sample_metrics.tsv (FP means), interim/viral_detections_long.tsv, interim/coverage/per_virus_coverage.tsv.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INT = ROOT / "analysis_708" / "interim"
COV = INT / "coverage"
PIPES = ["TaxProfiler", "TaxTriage", "mag-nf", "KrakenMetaReads-nf", "GOTTCHA2", "CLARK", "sourmash"]
PLATS = ["Illumina", "PacBio"]
THR = ["≥1", "≥10"]
THR_MIN = {"≥1": 1, "≥10": 10, "≥100": 100}
N_BOOT = 10000
rng = np.random.default_rng(20261012)

det = pd.read_csv(COV / "per_virus_detection_by_threshold.tsv", sep="\t")
cov = pd.read_csv(COV / "per_virus_coverage.tsv", sep="\t").set_index("taxid")
det["floor"] = det["taxid"].map(cov["pacbio_floor_applied"]).astype(bool)
vm = pd.read_csv(INT / "viral_sample_metrics.tsv", sep="\t")
vm = vm[vm.ground_truth == "non-phage virus (n=151)"]
fpmean = vm.groupby(["platform", "pipeline", "min_reads"])["FP"].mean()
n = len(det)
assert n == 151


def f1(recall, tp, fp):
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    return 2 * prec * recall / (prec + recall) if (prec + recall) > 0 else 0.0


# ------------------------------------------------------------------ 1. bootstrap over viruses
idx = rng.integers(0, n, size=(N_BOOT, n))
rows = []
for pipe in PIPES:
    for thr in THR:
        res = {}
        for plat in PLATS:
            x = det[f"{plat} | {pipe} | {thr}"].to_numpy() / 3.0           # per-virus detection frequency
            fp = fpmean[(plat, pipe, THR_MIN[thr])]
            rec_b = x[idx].mean(axis=1)
            f1_b = np.array([f1(r, r * n, fp) for r in rec_b])
            rec, tp = x.mean(), x.sum()
            res[plat] = (rec_b, f1_b)
            rows.append(dict(kind="level", pipeline=pipe, platform=plat, threshold=thr,
                             recall=rec, recall_lo=np.percentile(rec_b, 2.5), recall_hi=np.percentile(rec_b, 97.5),
                             f1=f1(rec, tp, fp), f1_lo=np.percentile(f1_b, 2.5), f1_hi=np.percentile(f1_b, 97.5)))
        d_rec = res["PacBio"][0] - res["Illumina"][0]
        d_f1 = res["PacBio"][1] - res["Illumina"][1]
        rec_pt = det[f"PacBio | {pipe} | {thr}"].mean() / 3 - det[f"Illumina | {pipe} | {thr}"].mean() / 3
        f1_pt = (f1(det[f"PacBio | {pipe} | {thr}"].mean() / 3, det[f"PacBio | {pipe} | {thr}"].sum() / 3,
                    fpmean[("PacBio", pipe, THR_MIN[thr])])
                 - f1(det[f"Illumina | {pipe} | {thr}"].mean() / 3, det[f"Illumina | {pipe} | {thr}"].sum() / 3,
                      fpmean[("Illumina", pipe, THR_MIN[thr])]))
        rows.append(dict(kind="difference (PacBio - Illumina)", pipeline=pipe, platform="", threshold=thr,
                         recall=rec_pt, recall_lo=np.percentile(d_rec, 2.5), recall_hi=np.percentile(d_rec, 97.5),
                         f1=f1_pt, f1_lo=np.percentile(d_f1, 2.5), f1_hi=np.percentile(d_f1, 97.5)))
boot = pd.DataFrame(rows)
boot.to_csv(COV / "bootstrap_ci.tsv", sep="\t", index=False)

# consistency with Table 3 (>=1)
agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t").set_index(["platform", "pipeline"])
chk = boot[(boot.kind == "level") & (boot.threshold == "≥1")]
diff = max(abs(r.recall - agg.loc[(r.platform, r.pipeline), "recall_mean"]) for r in chk.itertuples())
diff_f = max(abs(r.f1 - agg.loc[(r.platform, r.pipeline), "f1_mean"]) for r in chk.itertuples())
print(f"bootstrap point estimates vs Table 3: max |recall diff| = {diff:.2e}, max |F1 diff| = {diff_f:.3f}")

# ------------------------------------------------------------------ 2. metrics without floor-affected viruses
sub = det[~det["floor"]]
m = len(sub)
rows = []
for pipe in PIPES:
    for thr in THR:
        vals = {}
        for plat in PLATS:
            x = sub[f"{plat} | {pipe} | {thr}"].to_numpy() / 3.0
            rec, tp = x.mean(), x.sum()
            fp = fpmean[(plat, pipe, THR_MIN[thr])]
            prec = tp / (tp + fp) if tp + fp else 0.0
            vals[plat] = (rec, prec, f1(rec, tp, fp))
            rows.append(dict(pipeline=pipe, platform=plat, threshold=thr, n_viruses=m,
                             recall=rec, precision=prec, f1=vals[plat][2],
                             recall_all151=det[f"{plat} | {pipe} | {thr}"].mean() / 3))
        rows.append(dict(pipeline=pipe, platform="PacBio - Illumina", threshold=thr, n_viruses=m,
                         recall=vals["PacBio"][0] - vals["Illumina"][0], precision=vals["PacBio"][1] - vals["Illumina"][1],
                         f1=vals["PacBio"][2] - vals["Illumina"][2], recall_all151=np.nan))
pd.DataFrame(rows).to_csv(COV / "floor_excluded_metrics.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 3. false-positive profile
spec = importlib.util.spec_from_file_location("m14", ROOT / "analysis_708" / "scripts" / "14_viral_only_scoring.py")
m14 = importlib.util.module_from_spec(spec); sys.modules["m14"] = m14; spec.loader.exec_module(m14)
spec = importlib.util.spec_from_file_location("m15", ROOT / "analysis_708" / "scripts" / "15_genus_level_sensitivity.py")
m15 = importlib.util.module_from_spec(spec); sys.modules["m15"] = m15; spec.loader.exec_module(m15)
tax2genus, _, _ = m15.build_maps(m15.MAP_SOURCE)

gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
gt["taxid_num"] = pd.to_numeric(gt["taxid"], errors="coerce").astype("Int64")
npg = gt[gt.category.str.startswith("Virus (")]
ph = gt[gt.category == "Bacteriophage"]
truth_tax, truth_key = set(npg.taxid_num.dropna().astype(int)), set(npg.match_key)
ph_tax, ph_key = set(ph.taxid_num.dropna().astype(int)), set(ph.match_key)
truth_gen = {tax2genus.get(str(t), (None, None))[0] for t in truth_tax} - {None}

long = pd.read_csv(INT / "viral_detections_long.tsv", sep="\t", low_memory=False)
long = long[long["rank"].isin(m14.SPECIES_RANKS) & long["is_viral"].astype(str).eq("True")]
rows = []
for plat in PLATS:
    for pipe in PIPES:
        recs = []
        for sample, df in long[(long.platform == plat) & (long.pipeline == pipe)].groupby("sample"):
            uses_taxid = df["taxid"].dropna().astype(str).str.strip().ne("").any()
            s = df[df["reads"] > 0]
            if uses_taxid:
                s = s[~s.taxid_num.isin(truth_tax | ph_tax)].dropna(subset=["taxid_num"])
                s = s.groupby("taxid_num", as_index=False)["reads"].max()
                gen = s["taxid_num"].astype(int).astype(str).map(lambda t: tax2genus.get(t, (None, None))[0])
                s["related"] = gen.isin(truth_gen)
            else:
                s = s[~s.name_key.isin(truth_key | ph_key)]
                s = s.groupby("name_key", as_index=False)["reads"].max()
                s["related"] = np.nan
            recs.append(s.assign(sample=sample))
        fp = pd.concat(recs)
        per = fp.groupby("sample").size()
        rows.append(dict(platform=plat, pipeline=pipe, FP_nodes_mean=per.mean(),
                         pct_supported_by_1_read=100 * (fp.reads <= 1).mean(),
                         pct_supported_by_le3_reads=100 * (fp.reads <= 3).mean(),
                         pct_supported_by_ge10_reads=100 * (fp.reads >= 10).mean(),
                         median_reads=fp.reads.median(),
                         pct_in_genus_of_a_truth_virus=(100 * fp.related.mean() if fp.related.notna().any() else np.nan)))
fpp = pd.DataFrame(rows)
# GOTTCHA2 and sourmash report abundance-like values, not read counts, in the "reads" field: read-support statistics are
# comparable only for the Kraken2-based pipelines and CLARK
unit = fpp.pipeline.isin(["GOTTCHA2", "sourmash"])
fpp.loc[unit, ["pct_supported_by_1_read", "pct_supported_by_le3_reads", "pct_supported_by_ge10_reads", "median_reads"]] = np.nan
fpp["read_support_comparable"] = ~unit
fpp.to_csv(COV / "false_positive_profile.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 4. Kraken2 misses (Illumina) and database report
miss = det[det["Illumina | TaxProfiler | ≥1"] == 0].copy()
db_ids = {}
for line in open(ROOT / "Illumina/mag-nf/sample_0/Sample_708_illumina_0.kraken2_report.txt"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 6:
        db_ids[p[4].strip()] = int(p[1])
out = pd.DataFrame({
    "virus": miss["virus"], "taxid": miss["taxid"], "designed_relative_abundance_pct": miss["designed_relative_abundance_pct"],
    "simulated_read_pairs_sample0": miss["taxid"].map(lambda t: int(cov.loc[t, "illumina_pairs_s0"])),
    "taxid_present_in_kraken2_report": miss["taxid"].map(lambda t: str(t) in db_ids),
    "kraken2_clade_reads_sample0": miss["taxid"].map(lambda t: db_ids.get(str(t), np.nan)),
})
for p in PIPES:
    out[f"samples_detected_by_{p}"] = miss[f"Illumina | {p} | ≥1"].values
out.sort_values("designed_relative_abundance_pct", ascending=False).to_csv(COV / "kraken2_misses_illumina.tsv", sep="\t", index=False)
print("done")
print(pd.read_csv(COV / "false_positive_profile.tsv", sep="\t").round(1).to_string())
