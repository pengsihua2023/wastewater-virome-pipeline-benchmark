#!/usr/bin/env python3
"""Is viral recall bounded by classifier capability, or by how many reads each
virus actually contributed?

Recall over the full 151-member viral ground truth mixes two very different
situations: viruses that contributed enough reads to be findable, and viruses
that contributed essentially none. This script separates them.

For each virus, the reads it can be expected to contribute in a given library is

    expected_reads = (total viral reads in that library)
                     x (designed relative abundance of the virus / 4.00%)

where 4.00% is the viral block's share of the designed community. Detection is
then tabulated against that expectation, per pipeline and per platform.

CAVEAT for PacBio: the PBSIM3 driver applied a minimum-depth floor that raised
64.9% of viral members to roughly five reads irrespective of designed abundance
(see 10_simulation_config_audit.xlsx). Expected reads therefore *understate* the
signal actually present on PacBio, and the PacBio panel of this analysis should
be read as indicative only. The Illumina simulation applied no such floor.

Writes 13_detection_vs_depth.xlsx and interim/detection_vs_depth.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

VIRAL_BLOCK_SHARE = 4.00        # designed % of the community that is non-phage virus
BINS = [0, 0.5, 1, 2, 5, 20, np.inf]
LABELS = ["<0.5", "0.5-1", "1-2", "2-5", "5-20", ">20"]


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    v = gt[gt.category.str.startswith("Virus (")].copy()
    assert len(v) == 151, len(v)
    v["t"] = pd.to_numeric(v["taxid"], errors="coerce")

    comp = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
    viral_reads = comp.set_index(["platform", "sample"]).virus_reads.to_dict()

    det = pd.read_csv(INT / "viral_detection_detail.tsv", sep="\t")
    k2t = dict(zip(v.match_key, v.t))

    def as_taxid(k):
        try:
            return int(float(k))
        except (TypeError, ValueError):
            return k2t.get(str(k))

    det["t"] = det["key"].map(as_taxid)
    tp = det[det.status == "TP"]

    rows = []
    for (plat, samp), vr in viral_reads.items():
        exp = vr * v["relative_abundance_num"] / VIRAL_BLOCK_SHARE
        for pipe in sorted(det.pipeline.unique()):
            hit = set(tp[(tp.platform == plat) & (tp.pipeline == pipe) &
                         (tp["sample"] == samp)].t)
            d = pd.DataFrame({"t": v.t, "virus": v.std_name,
                              "abundance_pct": v.relative_abundance_num,
                              "expected_reads": exp.values})
            d["detected"] = d.t.isin(hit)
            d["platform"], d["pipeline"], d["sample"] = plat, pipe, samp
            rows.append(d)
    long = pd.concat(rows, ignore_index=True)
    long["depth_bin"] = pd.cut(long.expected_reads, BINS, right=False,
                               labels=LABELS).astype(str)
    long.to_csv(INT / "detection_vs_depth.tsv", sep="\t", index=False)

    # recall within each expected-depth bin
    binned = (long.groupby(["platform", "pipeline", "sample", "depth_bin"])
              .agg(n=("t", "size"), detected=("detected", "sum")).reset_index())
    binned["recall"] = binned.detected / binned.n
    agg = (binned.groupby(["platform", "pipeline", "depth_bin"])
           .agg(n_viruses=("n", "first"), recall=("recall", "mean")).reset_index())

    # recall restricted to viruses that could plausibly be found
    thr_rows = []
    for thr in [0, 1, 2, 5]:
        sub = long[long.expected_reads >= thr]
        g = (sub.groupby(["platform", "pipeline", "sample"])
             .agg(n=("t", "size"), detected=("detected", "sum")).reset_index())
        g["recall"] = g.detected / g.n
        gg = (g.groupby(["platform", "pipeline"])
              .agg(n_viruses=("n", "first"), recall=("recall", "mean")).reset_index())
        gg["min_expected_reads"] = thr
        thr_rows.append(gg)
    restricted = pd.concat(thr_rows, ignore_index=True)

    # viruses with real signal that were still missed by everything on a platform
    sig = long[long.expected_reads >= 5]
    missed_all = (sig.groupby(["platform", "t", "virus", "abundance_pct",
                               "expected_reads"]).detected.sum()
                  .reset_index(name="n_runs_detecting"))
    missed_all = missed_all[missed_all.n_runs_detecting == 0].sort_values(
        ["platform", "expected_reads"], ascending=[True, False])

    with pd.ExcelWriter(AN / "13_detection_vs_depth.xlsx") as xw:
        agg.to_excel(xw, sheet_name="recall_by_expected_depth", index=False)
        restricted.to_excel(xw, sheet_name="recall_above_threshold", index=False)
        missed_all.to_excel(xw, sheet_name="signal_present_but_missed", index=False)
        (long[long["sample"] == "sample_0"]
         [["platform", "pipeline", "virus", "abundance_pct", "expected_reads",
           "detected"]].to_excel(xw, sheet_name="per_virus_sample_0", index=False))

    ill = long[(long.platform == "Illumina") & (long["sample"] == "sample_0")]
    n_lt1 = (ill[ill.pipeline == "TaxProfiler"].expected_reads < 1).sum()
    n_lt2 = (ill[ill.pipeline == "TaxProfiler"].expected_reads < 2).sum()
    print(f"Illumina, per library: {int(ill.expected_reads.sum() / 1):,} expected viral reads "
          f"spread over 151 viruses")
    print(f"  viruses expected to contribute <1 read: {n_lt1}/151 ({n_lt1/151:.0%})")
    print(f"  viruses expected to contribute <2 reads: {n_lt2}/151\n")

    print("Recall restricted to viruses with at least N expected reads (mean of 3 samples):")
    piv = restricted.pivot_table(index=["platform", "pipeline"],
                                 columns="min_expected_reads", values="recall")
    print(piv.round(3).to_string())
    print("\nn viruses in each restricted set:")
    print(restricted.pivot_table(index="platform", columns="min_expected_reads",
                                 values="n_viruses", aggfunc="first").to_string())
    print("\nViruses with >=5 expected reads missed by ALL pipelines:")
    print(missed_all[["platform", "virus", "abundance_pct", "expected_reads"]]
          .round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
