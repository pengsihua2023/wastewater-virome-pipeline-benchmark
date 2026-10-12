#!/usr/bin/env python3
"""Verify that every headline number in the paper matches the interim tables.

Run last. Exits non-zero if any check fails.

Baseline updated 2026-09-28 for analysis version V4, after the simulated data was
regenerated to fix an abundance-unit bug (see guide-line-regenerat-data.md and
Verification-simulated-data-vs-design.md) and after a PacBio nf-core/mag
read-level-report duplicate-file bug was found and fixed (see 03_score.py and
14_viral_only_scoring.py). The V3 baseline this replaces remains available in
git history for anyone re-checking the pre-regeneration numbers.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
INT = ROOT / "analysis_708" / "interim"

CHECKS_RECALL = {
    ("Illumina", "mag-nf"): 0.9887, ("Illumina", "TaxProfiler"): 0.9873,
    ("Illumina", "TaxTriage"): 0.9873, ("Illumina", "KrakenMetaReads-nf"): 0.9779,
    ("Illumina", "CLARK"): 0.742, ("Illumina", "sourmash"): 0.5753,
    ("Illumina", "GOTTCHA2"): 0.5541, ("PacBio", "TaxProfiler"): 0.9831,
    ("PacBio", "TaxTriage"): 0.9798, ("PacBio", "CLARK"): 0.7354,
    ("PacBio", "mag-nf"): 0.9831, ("PacBio", "sourmash"): 0.3503,
    ("PacBio", "GOTTCHA2"): 0.5217, ("PacBio", "KrakenMetaReads-nf"): 0.0132,
}
CHECKS_OTHER = [
    ("Illumina", "GOTTCHA2", "precision", 0.781), ("Illumina", "GOTTCHA2", "f1", 0.6483),
    ("PacBio", "GOTTCHA2", "precision", 0.7971), ("PacBio", "mag-nf", "f1", 0.3132),
    ("PacBio", "TaxProfiler", "f1", 0.3132), ("PacBio", "GOTTCHA2", "f1", 0.6306),
]
CHECKS_FP = [
    ("Illumina", "GOTTCHA2", 110.0), ("Illumina", "KrakenMetaReads-nf", 5215.7),
    ("Illumina", "sourmash", 7174.7), ("Illumina", "TaxProfiler", 10236.0),
    ("Illumina", "mag-nf", 12038.0), ("Illumina", "CLARK", 17212.0),
    ("PacBio", "GOTTCHA2", 94.0), ("PacBio", "mag-nf", 3041.0),
    ("PacBio", "KrakenMetaReads-nf", 929.3),
]
# "non-phage viral recall" is the domain-stratified recall from the 708-member
# community-wide scoring (03_score.py -> 04_aggregate_errors_abundance.py's
# domain_recall.tsv), NOT the dedicated viral-only re-scoring in
# 14_viral_only_scoring.py that the paper's Table 3 reports. The two are
# related but not numerically identical, since community-wide scoring's FP
# accounting differs from the viral-only rules (Section 2.4 of the paper).
CHECKS_VIRAL_RECALL = [
    ("Illumina", "TaxProfiler", 0.9404), ("Illumina", "mag-nf", 0.947),
    ("Illumina", "CLARK", 0.5188), ("Illumina", "sourmash", 0.9073),
    ("PacBio", "TaxProfiler", 0.9272), ("PacBio", "TaxTriage", 0.9117),
    ("PacBio", "GOTTCHA2", 0.6402), ("PacBio", "mag-nf", 0.9272),
    ("PacBio", "sourmash", 0.0287),
]

fails = []


def chk(desc, ok, got=None):
    print(f"{'PASS' if ok else 'FAIL':6s} {desc}" + ("" if ok else f"  -> got {got}"))
    if not ok:
        fails.append(desc)


def main():
    agg = pd.read_csv(INT / "aggregate_metrics.tsv", sep="\t").set_index(["platform", "pipeline"])
    dom = pd.read_csv(INT / "domain_recall.tsv", sep="\t").set_index(["platform", "pipeline", "domain"])
    comp = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    miss = pd.read_csv(INT / "missed_summary.tsv", sep="\t")

    for k, v in CHECKS_RECALL.items():
        got = agg.loc[k, "recall_mean"]
        chk(f"{k[0]}/{k[1]} recall = {v}", abs(got - v) < 6e-4, round(got, 4))
    for plat, pipe, m, v in CHECKS_OTHER:
        got = agg.loc[(plat, pipe), m + "_mean"]
        chk(f"{plat}/{pipe} {m} = {v}", abs(got - v) < 6e-4, round(got, 4))
    for plat, pipe, v in CHECKS_FP:
        got = agg.loc[(plat, pipe), "FP_mean"]
        chk(f"{plat}/{pipe} FP = {v}", abs(got - v) < 1.0, round(got, 1))
    for plat, pipe, v in CHECKS_VIRAL_RECALL:
        got = dom.loc[(plat, pipe, "Virus (non-phage)"), "recall_mean"]
        chk(f"{plat}/{pipe} non-phage viral recall = {v}", abs(got - v) < 6e-4, round(got, 4))

    chk("max between-sample recall SD <= 0.0062", agg.recall_sd.max() <= 0.0062,
        round(agg.recall_sd.max(), 5))
    chk("recall range 0.013-0.989",
        abs(agg.recall_mean.min() - 0.013) < 1e-3 and abs(agg.recall_mean.max() - 0.989) < 1e-3)
    chk("precision range 0.010-0.797",
        abs(agg.precision_mean.min() - 0.010) < 1e-3 and abs(agg.precision_mean.max() - 0.797) < 1e-3)

    il = comp[comp.platform == "Illumina"].pct_virus
    pb = comp[comp.platform == "PacBio"].pct_virus
    chk("Illumina realised viral 11.30-11.52%",
        abs(il.min() - 11.30) < 5e-2 and abs(il.max() - 11.52) < 5e-2)
    chk("PacBio realised viral 11.67-11.78%",
        abs(pb.min() - 11.67) < 5e-2 and abs(pb.max() - 11.78) < 5e-2)

    chk("ground truth = 708 members", len(gt) == 708, len(gt))
    chk("split 490 / 67 / 151",
        (gt.category == "Bacteria/Archaea").sum() == 490
        and (gt.category == "Bacteriophage").sum() == 67
        and gt.category.str.startswith("Virus (").sum() == 151)
    chk("7 members missed by 36 of 42 runs", (miss.n_runs_missed == 36).sum() == 7,
        int((miss.n_runs_missed == 36).sum()))
    chk("no member missed by all 42 runs", (miss.n_runs_missed == 42).sum() == 0)

    print(f"\n{len(fails)} failures")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
