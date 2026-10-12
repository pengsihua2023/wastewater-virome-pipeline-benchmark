#!/usr/bin/env python3
"""Audit the simulation configuration against what the datasets actually contain.

Inputs now available:
  genome_locations.tsv   genome_ID -> FASTA path (taxid is the filename prefix,
                         the parent directory gives the block)
  distribution_{0,1,2}.txt  genome_ID -> abundance (sums to 1.0), one per sample

Both the CAMISIM Illumina runs and the PBSIM3 PacBio runs (simulate_pacbio.py)
were driven by these files.

The script:
  1. maps every genome_ID to its NCBI taxid and to the 708-member ground truth;
  2. checks the distribution files against the design document's abundances;
  3. estimates each member's genome length from the Illumina data, using the
     CAMISIM relation reads ~ copies x length, and validates it against known
     genome-size ranges;
  4. applies simulate_pacbio.py's own depth formula, including the
     --min_depth / --min_reads floor, and compares the predicted read counts
     with the PacBio data;
  5. quantifies how much of the read budget the floor consumes.

Writes 10_simulation_config_audit.xlsx.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AN = ROOT / "analysis_708"
INT = AN / "interim"

# parameters as used in run_pacbio_sample.sh
TOTAL_BASES = 5e9
MEAN_LENGTH = 15000
MIN_DEPTH = 1.0
MIN_READS = 5

ILLUMINA_TOTAL = {"sample_0": 16642897, "sample_1": 16643343, "sample_2": 16642601}
PACBIO_TOTAL = {"sample_0": 449134, "sample_1": 490009, "sample_2": 518608}


def main():
    # ------------------------------------------------ genome_locations mapping
    rows = []
    for line in open(ROOT / "genome_locations.tsv"):
        gid, path = line.rstrip("\n").split("\t")[:2]
        parts = path.split("/")
        block = parts[-2]
        m = re.match(r"^(\d+)_", parts[-1])
        rows.append(dict(genome_id=gid, block_dir=block,
                         taxid=int(m.group(1)) if m else np.nan,
                         fasta=parts[-1]))
    loc = pd.DataFrame(rows)
    print(f"genome_locations.tsv: {len(loc)} genomes, "
          f"{loc.taxid.notna().sum()} with a parseable taxid")
    print(loc.block_dir.value_counts().to_string())

    # ------------------------------------------------------------ ground truth
    gt = pd.read_csv(INT / "ground_truth_raw.tsv", sep="\t")
    gt["taxid"] = pd.to_numeric(gt["taxid"], errors="coerce")
    gt["blk"] = np.where(gt.category == "Bacteria/Archaea", "Bacteria_Archaea",
                np.where(gt.category == "Bacteriophage", "Bacteriophage", "Viruses"))
    m = loc.merge(gt[["taxid", "std_name", "blk", "relative_abundance_num",
                      "design_reads_num", "category"]], on="taxid", how="left")
    print(f"\nmapped to ground truth: {m.std_name.notna().sum()} / {len(m)}")
    mismatch = m[(m.blk.notna()) & (m.block_dir != m.blk)]
    print(f"block directory disagrees with ground-truth category: {len(mismatch)}")

    # --------------------------------------------------------- distributions
    for i in range(3):
        d = pd.read_csv(ROOT / f"distribution_{i}.txt", sep="\t", header=None,
                        names=["genome_id", f"abund_{i}"])
        m = m.merge(d, on="genome_id", how="left")
    m["abund_mean"] = m[[f"abund_{i}" for i in range(3)]].mean(axis=1)

    # does distribution_0 reproduce the design document's abundance column?
    cmp = m.dropna(subset=["relative_abundance_num"])
    delta = (100 * cmp["abund_0"] - cmp["relative_abundance_num"]).abs()
    print(f"\ndistribution_0 vs design document abundance column:")
    print(f"  max |difference| = {delta.max():.8f} percentage points "
          f"({'identical' if delta.max() < 1e-4 else 'DIFFERENT'})")

    # --------------------- genome length estimated from the Illumina datasets
    # CAMISIM: reads_i  ~  copies_i x length_i  ->  length_i  ~  reads_i / abundance_i
    det = pd.read_csv(INT / "detections_long.tsv", sep="\t",
                      dtype={"taxid": str}, low_memory=False)
    det["t"] = pd.to_numeric(det["taxid"], errors="coerce")
    ill = det[(det.platform == "Illumina") & (det.parser == "taxprofiler_kraken2")]
    pb = det[(det.platform == "PacBio") & (det.parser == "taxprofiler_kraken2")]
    for samp in ILLUMINA_TOTAL:
        m[f"ill_{samp}"] = m.taxid.map(ill[ill["sample"] == samp].groupby("t").reads.max())
        m[f"pb_{samp}"] = m.taxid.map(pb[pb["sample"] == samp].groupby("t").reads.max())
    m["ill_reads"] = m[[f"ill_{s}" for s in ILLUMINA_TOTAL]].mean(axis=1).fillna(0)
    m["pb_reads"] = m[[f"pb_{s}" for s in PACBIO_TOTAL]].mean(axis=1).fillna(0)

    ok = m[(m.ill_reads > 0) & (m.abund_mean > 0)].copy()
    ok["len_proxy"] = ok["ill_reads"] / ok["abund_mean"]
    # calibrate so that the bacterial/archaeal median equals 4.0 Mb
    scale = 4.0e6 / ok.loc[ok.block_dir == "Bacteria_Archaea", "len_proxy"].median()
    m["est_genome_bp"] = (m["ill_reads"] / m["abund_mean"]) * scale
    # members with zero Illumina reads yield no estimate; fill with the median of
    # the members of their own block that were detected, so the floor
    # calculation below is defined for every genome
    m.loc[m["ill_reads"] <= 0, "est_genome_bp"] = np.nan
    blockmed = m.groupby("block_dir").est_genome_bp.median()
    m["est_genome_bp_source"] = np.where(m.est_genome_bp.notna(),
                                         "estimated from Illumina", "block median")
    m["est_genome_bp"] = m["est_genome_bp"].fillna(m.block_dir.map(blockmed))
    est = blockmed
    print("\nGenome length estimated from the Illumina data, median over members "
          "with >0 Illumina reads (calibrated: bacterial median = 4.0 Mb):")
    print(est.round(0).to_string())
    print("  n with a direct estimate per block:")
    print(m[m.est_genome_bp_source == "estimated from Illumina"]
          .block_dir.value_counts().to_string())
    print("  -> compare with expected: Bacteria ~4 Mb, ssRNA phage ~4 kb, "
          "small RNA viruses ~7 kb")

    # ------------------- apply simulate_pacbio.py's own formula to each genome
    L = m["est_genome_bp"]
    depth_from_abundance = TOTAL_BASES * m["abund_mean"] / L
    floor_reads = (MIN_READS * MEAN_LENGTH) / L
    effective_floor = np.maximum(MIN_DEPTH, floor_reads)
    m["depth_requested"] = depth_from_abundance
    m["depth_floor"] = effective_floor
    m["depth_applied"] = np.maximum(depth_from_abundance, effective_floor)
    m["floored"] = depth_from_abundance < effective_floor
    # reads a genome yields: a read cannot exceed the template length
    eff_read_len = np.minimum(MEAN_LENGTH, L)
    m["pb_reads_predicted"] = m["depth_applied"] * L / eff_read_len
    m["pb_ratio_obs_over_pred"] = m["pb_reads"] / m["pb_reads_predicted"]

    print("\nsimulate_pacbio.py depth floor (--min_depth 1.0, --min_reads 5):")
    fl = m.groupby("block_dir").floored.agg(["sum", "size", "mean"])
    fl.columns = ["n_floored", "n_total", "fraction_floored"]
    print(fl.round(3).to_string())
    budget = (m["depth_applied"] * L).groupby(m.block_dir).sum()
    req = (depth_from_abundance * L).groupby(m.block_dir).sum()
    print(f"\nbases implied by the floored depths: {budget.sum()/1e9:.2f} Gbp "
          f"vs --total_bases 5.00 Gbp  (over-run {budget.sum()/TOTAL_BASES:.2f}x)")
    print((pd.DataFrame({"requested_Gbp": req / 1e9, "after_floor_Gbp": budget / 1e9,
                         "inflation": budget / req}).round(3)).to_string())

    print("\nObserved / predicted PacBio reads, by block (median):")
    print(m.groupby("block_dir").pb_ratio_obs_over_pred.median().round(5).to_string())

    with pd.ExcelWriter(AN / "10_simulation_config_audit.xlsx") as xw:
        cols = ["genome_id", "taxid", "std_name", "block_dir", "category", "fasta",
                "abund_0", "abund_1", "abund_2", "abund_mean",
                "relative_abundance_num", "design_reads_num", "est_genome_bp",
                "depth_requested", "depth_floor", "depth_applied", "floored",
                "est_genome_bp_source", "ill_reads", "pb_reads", "pb_reads_predicted",
                "pb_ratio_obs_over_pred"]
        m[cols].to_excel(xw, sheet_name="per_genome", index=False)
        fl.reset_index().to_excel(xw, sheet_name="floor_summary", index=False)
        m.loc[m.floored, cols].to_excel(xw, sheet_name="floored_genomes", index=False)
        m.loc[m.block_dir != "Bacteria_Archaea", cols].to_excel(
            xw, sheet_name="viral_genomes", index=False)

    m.to_csv(INT / "simulation_config_audit.tsv", sep="\t", index=False)
    print("\nwrote 10_simulation_config_audit.xlsx")
    return 0


if __name__ == "__main__":
    sys.exit(main())
