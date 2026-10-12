#!/usr/bin/env python3
"""Per-member comparison of designed read counts against realised read counts.

The design document states per-member read counts explicitly ("Simulation design
reads", equal to the actual Kraken2 clade reads of sample HTP-2024-01-15 for all
automatically selected members). This script compares those numbers directly
with what the simulated datasets actually delivered, avoiding any inference
about abundance units.

It also discriminates which of the document's two mutually inconsistent
specifications was passed to CAMISIM:
  (a) the "Simulation design reads" column  -> viral share 93.4% of reads
  (b) the "Relative abundance in simulation (%)" column (88/8/4) -> 12%
by predicting the viral read share each would produce once weighted by genome
length, and comparing with the measured value.

Writes 09_designed_vs_observed_reads.xlsx.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"

# order-of-magnitude genome sizes, for the discrimination test only
APPROX_GENOME_BP = {"Bacteria/Archaea": 4.0e6, "Bacteriophage": 4.0e3, "Virus": 7.0e3}

STREAM = {"Illumina": "taxprofiler_kraken2", "PacBio": "taxprofiler_kraken2"}
# total reads per sample, read live from interim/realised_composition.tsv
# (05_realised_composition.py) rather than hand-copied, so this never goes
# stale when the underlying data is regenerated.
_rc = pd.read_csv(INT / "realised_composition.tsv", sep="\t")
TOTALS = {(r.platform, r.sample): int(r.total_reads) for r in _rc.itertuples()}


def main():
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gt["blk"] = np.where(gt.category == "Bacteria/Archaea", "Bacteria/Archaea",
                np.where(gt.category == "Bacteriophage", "Bacteriophage", "Virus"))
    gt["t"] = pd.to_numeric(gt["taxid"], errors="coerce")

    det = pd.read_csv(INT / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    det["t"] = pd.to_numeric(det["taxid"], errors="coerce")

    per_member = gt[["blk", "category", "raw_name", "std_name", "taxid", "t",
                     "design_reads_num", "relative_abundance_num"]].copy()
    per_member = per_member.rename(columns={"design_reads_num": "designed_reads"})

    design_total = per_member["designed_reads"].sum()
    for plat, parser in STREAM.items():
        for samp in ["sample_0", "sample_1", "sample_2"]:
            s = det[(det.platform == plat) & (det.parser == parser) &
                    (det["sample"] == samp)]
            obs = s.groupby("t").reads.max()
            col = f"observed_{plat}_{samp[-1]}"
            per_member[col] = per_member["t"].map(obs).fillna(0.0)
        cols = [f"observed_{plat}_{i}" for i in range(3)]
        per_member[f"observed_{plat}_mean"] = per_member[cols].mean(axis=1)
        # what the member should have received if the design read counts had been
        # honoured proportionally within the realised library size
        scale = np.mean([TOTALS[(plat, f"sample_{i}")] for i in range(3)]) / design_total
        per_member[f"expected_{plat}_if_design_honoured"] = per_member["designed_reads"] * scale
        per_member[f"shortfall_{plat}_x"] = (
            per_member[f"expected_{plat}_if_design_honoured"] /
            per_member[f"observed_{plat}_mean"].replace(0, np.nan))

    # ---------------------------------------------------------- block summary
    blk = per_member.groupby("blk").agg(
        n_members=("std_name", "size"),
        designed_reads=("designed_reads", "sum"),
        designed_read_share_pct=("designed_reads", lambda x: 100 * x.sum() / design_total),
        design_abundance_pct=("relative_abundance_num", "sum"),
        observed_Illumina=("observed_Illumina_mean", "sum"),
        observed_PacBio=("observed_PacBio_mean", "sum"))
    for plat in ["Illumina", "PacBio"]:
        tot = np.mean([TOTALS[(plat, f"sample_{i}")] for i in range(3)])
        blk[f"observed_read_share_pct_{plat}"] = 100 * blk[f"observed_{plat}"] / tot
        blk[f"delivered_per_designed_read_{plat}"] = blk[f"observed_{plat}"] / blk["designed_reads"]
    for plat in ["Illumina", "PacBio"]:
        base = blk.loc["Bacteria/Archaea", f"delivered_per_designed_read_{plat}"]
        blk[f"delivery_vs_bacteria_{plat}"] = blk[f"delivered_per_designed_read_{plat}"] / base

    # ------------------------------------------- which column drove CAMISIM?
    disc = []
    for label, col in [("Simulation design reads (viral 93.4% of reads)", "designed_reads"),
                       ("Relative abundance 88/8/4 (viral 12%)", "relative_abundance_num")]:
        w = per_member.groupby("blk")[col].sum()
        pred = {k: w[k] * APPROX_GENOME_BP[k] for k in APPROX_GENOME_BP}
        share = 100 * (pred["Bacteriophage"] + pred["Virus"]) / sum(pred.values())
        disc.append(dict(column_hypothesised_as_CAMISIM_input=label,
                         predicted_viral_read_share_pct=share))
    disc.append(dict(column_hypothesised_as_CAMISIM_input="MEASURED (Illumina)",
                     predicted_viral_read_share_pct=100 * blk.loc[["Bacteriophage", "Virus"],
                     "observed_Illumina"].sum() / np.mean([TOTALS[("Illumina", f"sample_{i}")]
                                                           for i in range(3)])))
    disc = pd.DataFrame(disc)

    with pd.ExcelWriter(AN / "09_designed_vs_observed_reads.xlsx") as xw:
        per_member.sort_values("designed_reads", ascending=False).to_excel(
            xw, sheet_name="per_member", index=False)
        blk.reset_index().to_excel(xw, sheet_name="by_block", index=False)
        disc.to_excel(xw, sheet_name="which_column_drove_camisim", index=False)
        per_member[per_member.blk == "Virus"].sort_values(
            "designed_reads", ascending=False).to_excel(
            xw, sheet_name="viruses_only", index=False)

    per_member.to_csv(INT / "designed_vs_observed_reads.tsv", sep="\t", index=False)

    print(blk[["n_members", "designed_reads", "designed_read_share_pct",
               "design_abundance_pct", "observed_read_share_pct_Illumina",
               "delivery_vs_bacteria_Illumina",
               "delivery_vs_bacteria_PacBio"]].round(6).to_string())
    print()
    print(disc.round(4).to_string(index=False))
    print()
    v = per_member[per_member.blk == "Virus"].nlargest(12, "designed_reads")
    print(v[["std_name", "designed_reads", "expected_Illumina_if_design_honoured",
             "observed_Illumina_mean", "shortfall_Illumina_x"]].round(1).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
