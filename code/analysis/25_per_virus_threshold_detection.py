#!/usr/bin/env python3
"""Per-virus detection of the 151 non-phage ground-truth viruses at the three read-count thresholds.

Re-applies the rules of 13_viral_only_scoring.py (species-or-below nodes, exact taxid match, sourmash matched by
normalised name, node kept when reads > 0 / 9 / 99, i.e. >=1 / >=10 / >=100 clade reads) to the stored
interim/viral_detections_long.tsv, and checks the result against interim/viral_sample_metrics.tsv.

Output: interim/coverage/per_virus_detection_by_threshold.tsv -- one row per virus; for every
platform x pipeline x threshold the number of the three samples (0-3) in which the virus was detected,
plus the number of pipelines detecting it in at least one sample.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
INT = ROOT / "analysis_708" / "interim"
spec = importlib.util.spec_from_file_location("m14", REPO / "code" / "analysis" / "13_viral_only_scoring.py")
m14 = importlib.util.module_from_spec(spec)
sys.modules["m14"] = m14
spec.loader.exec_module(m14)

PIPES = ["TaxProfiler", "TaxTriage", "mag-nf", "KrakenMetaReads-nf", "GOTTCHA2", "CLARK", "sourmash"]
THR = {"≥1": 0, "≥10": 9, "≥100": 99}

gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
gt["taxid_num"] = pd.to_numeric(gt["taxid"], errors="coerce").astype("Int64")
nonphage = gt[gt.category.str.startswith("Virus (")].copy()
assert len(nonphage) == 151
cov = pd.read_csv(INT / "coverage" / "per_virus_coverage.tsv", sep="\t").set_index("taxid")

long = pd.read_csv(INT / "viral_detections_long.tsv", sep="\t", low_memory=False)
long = long[long["rank"].isin(m14.SPECIES_RANKS)]

det = {}          # (platform, pipeline, thr_label) -> Series(index=taxid) of sample counts
per_run_tp = {}   # for validation
for platform in ["Illumina", "PacBio"]:
    for pipe in PIPES:
        runs = long[(long.platform == platform) & (long.pipeline == pipe)]
        counts = {lab: pd.Series(0, index=nonphage.taxid_num.astype(int)) for lab in THR}
        for sample, df in runs.groupby("sample"):
            uses_taxid = df["taxid"].dropna().astype(str).str.strip().ne("").any()
            for lab, thr in THR.items():
                s = df[df["reads"] > thr]
                if uses_taxid:
                    hit = set(s["taxid_num"].dropna().astype(int))
                    flag = nonphage.taxid_num.astype(int).isin(hit)
                else:
                    hit = set(s["name_key"]) - {""}
                    flag = nonphage.match_key.isin(hit)
                flag.index = nonphage.taxid_num.astype(int)
                counts[lab] = counts[lab] + flag.astype(int).values
                per_run_tp[(platform, pipe, sample, thr + 1)] = int(flag.sum())
        for lab in THR:
            det[(platform, pipe, lab)] = counts[lab]

# ---- validation against the stored per-run metrics (TP counts)
vm = pd.read_csv(INT / "viral_sample_metrics.tsv", sep="\t")
vm = vm[vm.ground_truth == "non-phage virus (n=151)"]
bad = 0
for r in vm.itertuples():
    mine = per_run_tp[(r.platform, r.pipeline, r.sample, r.min_reads)]
    if mine != r.TP:
        bad += 1
        print("MISMATCH", r.platform, r.pipeline, r.sample, r.min_reads, mine, r.TP)
assert bad == 0, bad
print("validated: per-virus detections reproduce all", len(vm), "TP counts in viral_sample_metrics.tsv")

out = nonphage[["taxid_num", "std_name", "category", "relative_abundance_num"]].copy()
out.columns = ["taxid", "virus", "category", "designed_relative_abundance_pct"]
out["taxid"] = out["taxid"].astype(int)
out["sentinel"] = out["category"].str.startswith("Virus (manual")
out = out.set_index("taxid")
out["illumina_reads_mean"] = cov["illumina_reads_mean"]
out["pacbio_reads_mean"] = cov["pacbio_reads_mean"]
for platform in ["Illumina", "PacBio"]:
    for lab in THR:
        anyp = sum((det[(platform, p, lab)] > 0).astype(int) for p in PIPES)
        out[f"{platform}: pipelines detecting (of 7), {lab}"] = anyp
for platform in ["Illumina", "PacBio"]:
    for pipe in PIPES:
        for lab in THR:
            out[f"{platform} | {pipe} | {lab}"] = det[(platform, pipe, lab)]
out = out.reset_index().sort_values(["sentinel", "designed_relative_abundance_pct"], ascending=[False, False])
out.to_csv(INT / "coverage" / "per_virus_detection_by_threshold.tsv", sep="\t", index=False)
print("wrote", INT / "coverage" / "per_virus_detection_by_threshold.tsv", out.shape)

# summary used in the paper text
for platform in ["Illumina", "PacBio"]:
    for lab in THR:
        c = out[f"{platform}: pipelines detecting (of 7), {lab}"]
        s = out[out.sentinel][f"{platform}: pipelines detecting (of 7), {lab}"]
        print(f"{platform} {lab}: viruses detected by >=1 pipeline {int((c>0).sum())}/151; by >=3 pipelines {int((c>=3).sum())}; "
              f"sentinels detected by >=1 pipeline {int((s>0).sum())}/13")
