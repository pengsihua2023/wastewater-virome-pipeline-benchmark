#!/usr/bin/env python3
"""Complete simulated community (708 members) with designed abundance and simulated read counts.

One row per community member: 151 non-phage viruses (138 automatically selected, 13 manually added sentinels),
67 bacteriophages and 490 bacteria/archaea. Read counts are taken from the read headers of the regenerated
FASTQ files (interim/coverage/raw_*): Illumina read pairs (R1 count) and PacBio HiFi reads, per sample, together
with the mean number of bases and mean depth.

Output: interim/coverage/full_community_708.tsv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COV = ROOT / "analysis_708" / "interim" / "coverage"
REGEN = ROOT / "regen"
S = [0, 1, 2]

gl = pd.read_csv(REGEN / "genome_lengths.tsv", sep="\t")
gt = pd.read_csv(ROOT / "analysis_708" / "interim" / "ground_truth_raw.tsv", sep="\t")
gt["taxid_num"] = pd.to_numeric(gt["taxid"], errors="coerce").astype("Int64")
assert gt.taxid_num.is_unique
info = gt.set_index("taxid_num")[["std_name", "ncbi_rank", "category"]]
gl = gl.join(info, on="taxid")
assert gl.std_name.notna().all()

GROUP = {"Viruses": "non-phage virus", "Bacteriophage": "bacteriophage", "Bacteria_Archaea": "bacterium / archaeon"}
gl["group"] = gl["dom"].map(GROUP)
gl.loc[gl.category.fillna("").str.startswith("Virus (manual"), "group"] = "non-phage virus (manually added sentinel)"

acc2gid = {}
for gid, path in zip(gl["gid"], gl["path"]):
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                acc2gid[line[1:].split()[0]] = gid

out = gl[["gid", "taxid", "std_name", "group", "ncbi_rank", "length_bp", "n_contigs"]].rename(
    columns={"std_name": "species", "length_bp": "genome_length_bp"}).set_index("gid")
plan = pd.read_csv(REGEN / "plan_pacbio_s0.tsv", sep="\t").set_index("gid")
out["designed_share_pct"] = 100 * plan["ab"]
pi = pd.read_csv(REGEN / "plan_illumina.tsv", sep="\t").set_index("gid")
out["designed_illumina_pairs_s0"] = pi["target_pairs"]

ib, pb = [], []
for s in S:
    parts = []
    for r in ("01", "02"):
        df = pd.read_csv(COV / f"raw_illumina_s{s}_r{r}_peracc.tsv", sep="\t", header=None, names=["acc", "reads", "bases"])
        df["gid"] = df["acc"].map(acc2gid)
        assert df.gid.notna().all()
        parts.append(df.groupby("gid")[["reads", "bases"]].sum())
    out[f"illumina_pairs_s{s}"] = parts[0]["reads"].reindex(out.index).fillna(0).astype(int)
    ib.append((parts[0]["bases"] + parts[1]["bases"]).reindex(out.index).fillna(0))
    p = pd.read_csv(COV / f"raw_pacbio_s{s}_pergenome.tsv", sep="\t", header=None, names=["gid", "reads", "bases"]).set_index("gid")
    out[f"pacbio_reads_s{s}"] = p["reads"].reindex(out.index).fillna(0).astype(int)
    pb.append(p["bases"].reindex(out.index).fillna(0))
out["illumina_bases_mean"] = sum(ib) / 3
out["pacbio_bases_mean"] = sum(pb) / 3
out["illumina_depth_x"] = out["illumina_bases_mean"] / out["genome_length_bp"]
out["pacbio_depth_x"] = out["pacbio_bases_mean"] / out["genome_length_bp"]

# ---- reads in the real source sample HTP-2024-01-15 and the design reads used for the simulation
src = gt.set_index("taxid_num")
out = out.reset_index()
out["source_reads_direct"] = out["taxid"].map(src["actual_reads_direct"])
out["source_reads_clade"] = out["taxid"].map(src["actual_reads_clade"])
out["design_reads"] = out["taxid"].map(src["design_reads_num"]).astype(int)
cat = out["taxid"].map(src["category"])
basis = pd.Series("equals source-sample clade reads", index=out.index)
basis[cat == "Virus (manual, no-rank)"] = "equals source-sample direct reads"
basis[cat == "Virus (manual, low-read)"] = "manually assigned (below the 20-read filter in the source sample)"
basis[cat == "Virus (manual, not detected)"] = "manually assigned (not detected in the source sample)"
out["design_basis"] = basis
nd = cat == "Virus (manual, not detected)"
out.loc[nd, ["source_reads_direct", "source_reads_clade"]] = 0
# design reads are weights within each domain; domains are fixed at 4% (viruses), 8% (phages), 88% (bacteria/archaea)
frac = {"non-phage virus": 4.0, "non-phage virus (manually added sentinel)": 4.0, "bacteriophage": 8.0, "bacterium / archaeon": 88.0}
dom = out["group"].map(lambda g: "virus" if g.startswith("non-phage") else g)
tot = out.groupby(dom)["design_reads"].transform("sum")
check = out["group"].map(frac) * out["design_reads"] / tot
assert (check - out["designed_share_pct"]).abs().max() < 1e-3, (check - out["designed_share_pct"]).abs().max()
print("designed_share_pct = domain fraction x design_reads / sum of design_reads within the domain (max abs diff %.2e)" %
      (check - out["designed_share_pct"]).abs().max())

out.insert(out.columns.get_loc("group") + 1, "sentinel_virus", out["group"].str.contains("sentinel"))
order = {"non-phage virus (manually added sentinel)": 0, "non-phage virus": 1, "bacteriophage": 2, "bacterium / archaeon": 3}
out["_o"] = out["group"].map(order)
out = out.sort_values(["_o", "designed_share_pct"], ascending=[True, False]).drop(columns="_o")
out.to_csv(COV / "full_community_708.tsv", sep="\t", index=False)
print(out.shape, out.group.value_counts().to_dict())
tot = out[[f"illumina_pairs_s{s}" for s in S]].sum()
print("Illumina read pairs per sample:", tot.tolist(), "| PacBio reads per sample:", out[[f"pacbio_reads_s{s}" for s in S]].sum().tolist())
print("members with zero reads (Illumina, PacBio):", int((out[[f'illumina_pairs_s{s}' for s in S]].sum(axis=1) == 0).sum()),
      int((out[[f'pacbio_reads_s{s}' for s in S]].sum(axis=1) == 0).sum()))
