#!/usr/bin/env python3
"""Per-class and per-virus read/base coverage of the six simulated datasets.

Inputs (all produced from the regenerated FASTQ files in regen/):
  interim/coverage/raw_illumina_s{0,1,2}_r{01,02}_peracc.tsv   accession, reads, bases
  interim/coverage/raw_pacbio_s{0,1,2}_pergenome.tsv            genome id, reads, bases
  regen/plan_pacbio_s{0,1,2}.tsv, regen/plan_illumina*.tsv      designed targets
  regen/genome_lengths.tsv, interim/ground_truth_raw.tsv

Outputs (interim/coverage/):
  class_summary.tsv, per_virus_coverage.tsv, sentinel_coverage.tsv,
  pacbio_read_qc.tsv, floor_summary.tsv
"""
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COV = ROOT / "analysis_708" / "interim" / "coverage"
REGEN = ROOT / "regen"
SAMPLES = [0, 1, 2]

gl = pd.read_csv(REGEN / "genome_lengths.tsv", sep="\t")
gt = pd.read_csv(ROOT / "analysis_708" / "interim" / "ground_truth_raw.tsv", sep="\t")
gt["taxid_num"] = gt["taxid_num"].astype(int)
cat = gt.set_index("taxid_num")["category"].to_dict()
name = gt.set_index("taxid_num")["std_name"].to_dict()
gl["category"] = gl["taxid"].map(cat)
gl["name"] = gl["taxid"].map(name)
gl["klass"] = gl["dom"].map({"Viruses": "non-phage viruses (151 targets)",
                             "Bacteriophage": "bacteriophages (67)",
                             "Bacteria_Archaea": "bacteria and archaea (490)"})
gl["sentinel"] = gl["category"].fillna("").str.startswith("Virus (manual")

# ---- accession -> genome id (Illumina read headers carry the source accession)
acc2gid = {}
for gid, path in zip(gl["gid"], gl["path"]):
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                acc2gid[line[1:].split()[0]] = gid

ill = {}
for s in SAMPLES:
    tot = {}
    for r in ("01", "02"):
        df = pd.read_csv(COV / f"raw_illumina_s{s}_r{r}_peracc.tsv", sep="\t", header=None,
                         names=["acc", "reads", "bases"])
        df["gid"] = df["acc"].map(acc2gid)
        assert df["gid"].notna().all(), df[df["gid"].isna()].head()
        tot[r] = df.groupby("gid")[["reads", "bases"]].sum()
    pairs = tot["01"]["reads"]
    assert (pairs == tot["02"]["reads"]).all()
    ill[s] = pd.DataFrame({"reads": pairs, "bases": tot["01"]["bases"] + tot["02"]["bases"]})

pac = {}
for s in SAMPLES:
    df = pd.read_csv(COV / f"raw_pacbio_s{s}_pergenome.tsv", sep="\t", header=None,
                     names=["gid", "reads", "bases"]).set_index("gid")
    pac[s] = df

# ---- class summary
rows = []
for plat, data in (("Illumina (read pairs; bases R1+R2)", ill), ("PacBio HiFi", pac)):
    for s in SAMPLES:
        d = data[s].join(gl.set_index("gid")["klass"])
        tr, tb = d["reads"].sum(), d["bases"].sum()
        for k, g in d.groupby("klass"):
            rows.append(dict(platform=plat, sample=s, klass=k, reads=int(g["reads"].sum()),
                             bases=int(g["bases"].sum()), read_share_pct=100 * g["reads"].sum() / tr,
                             base_share_pct=100 * g["bases"].sum() / tb))
cls = pd.DataFrame(rows)
cls.to_csv(COV / "class_summary.tsv", sep="\t", index=False)

# ---- per-virus coverage
v = gl[gl["dom"] == "Viruses"].copy().set_index("gid")
plan = {s: pd.read_csv(REGEN / f"plan_pacbio_s{s}.tsv", sep="\t").set_index("gid") for s in SAMPLES}
pi = pd.read_csv(REGEN / "plan_illumina.tsv", sep="\t").set_index("gid")
out = v[["taxid", "name", "category", "sentinel", "length_bp", "n_contigs"]].copy()
out["designed_share_pct"] = [100 * plan[0].loc[g, "ab"] for g in out.index]
out["designed_illumina_pairs_s0"] = pi.loc[out.index, "target_pairs"]
for s in SAMPLES:
    T = plan[s]["target_zmw"].sum()
    out[f"pacbio_expected_reads_s{s}"] = [plan[s].loc[g, "ab"] * T for g in out.index]
    out[f"pacbio_target_zmw_s{s}"] = [plan[s].loc[g, "target_zmw"] for g in out.index]
    out[f"illumina_pairs_s{s}"] = ill[s]["reads"].reindex(out.index).fillna(0).astype(int)
    out[f"pacbio_reads_s{s}"] = pac[s]["reads"].reindex(out.index).fillna(0).astype(int)
    out[f"illumina_bases_s{s}"] = ill[s]["bases"].reindex(out.index).fillna(0).astype(int)
    out[f"pacbio_bases_s{s}"] = pac[s]["bases"].reindex(out.index).fillna(0).astype(int)
for plat, unit in (("illumina", "pairs"), ("pacbio", "reads")):
    col = "illumina_pairs" if plat == "illumina" else "pacbio_reads"
    out[f"{plat}_reads_mean"] = out[[f"{col}_s{s}" for s in SAMPLES]].mean(axis=1)
    out[f"{plat}_reads_min"] = out[[f"{col}_s{s}" for s in SAMPLES]].min(axis=1)
    out[f"{plat}_bases_mean"] = out[[f"{plat}_bases_s{s}" for s in SAMPLES]].mean(axis=1)
    out[f"{plat}_depth_x"] = out[f"{plat}_bases_mean"] / out["length_bp"]
out["pacbio_expected_mean"] = out[[f"pacbio_expected_reads_s{s}" for s in SAMPLES]].mean(axis=1)
out["pacbio_floor_applied"] = out["pacbio_expected_mean"] < 1
out = out.reset_index().sort_values("designed_share_pct", ascending=False)
out.to_csv(COV / "per_virus_coverage.tsv", sep="\t", index=False)
out[out["sentinel"]].to_csv(COV / "sentinel_coverage.tsv", sep="\t", index=False)

# ---- floor summary
fl = []
for s in SAMPLES:
    e = out[f"pacbio_expected_reads_s{s}"]
    t = out[f"pacbio_target_zmw_s{s}"]
    fl.append(dict(sample=s, viruses_target_eq_1=int((t <= 1).sum()), viruses_expected_lt_1=int((e < 1).sum()),
                   min_expected=float(e.min()),
                   min_observed_pacbio=int(out[f"pacbio_reads_s{s}"].min()),
                   min_observed_illumina=int(out[f"illumina_pairs_s{s}"].min()),
                   illumina_min_designed_pairs=int(out["designed_illumina_pairs_s0"].min())))
pd.DataFrame(fl).to_csv(COV / "floor_summary.tsv", sep="\t", index=False)

# ---- PacBio read length / accuracy
qc = []
for s in SAMPLES:
    h = pd.read_csv(COV / f"raw_pacbio_s{s}_lenhist.tsv", sep="\t", header=None, names=["bin", "n"]).sort_values("bin")
    cum = h["n"].cumsum() / h["n"].sum()
    q = {p: int(h.loc[cum >= p / 100, "bin"].iloc[0]) for p in (5, 25, 50, 75, 95)}
    n, tot = map(int, open(COV / f"raw_pacbio_s{s}_total.tsv").read().split())
    # accuracy from the first 100,000 reads
    accs = []
    with gzip.open(REGEN / "pacbio" / f"pacbio_hifi_reads_sample_{s}.fastq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                perr = np.mean([10 ** (-(ord(c) - 33) / 10) for c in line.strip()])
                accs.append(1 - perr)
                if len(accs) >= 100000:
                    break
    accs = np.array(accs)
    qc.append(dict(sample=s, reads=n, total_bases=tot, mean_length=tot / n,
                   len_p5_bin=q[5], len_p25_bin=q[25], len_median_bin=q[50], len_p75_bin=q[75], len_p95_bin=q[95],
                   accuracy_mean_first100k=accs.mean(), accuracy_median=np.median(accs),
                   frac_ge_Q20=float((accs >= 0.99).mean()), frac_ge_Q30=float((accs >= 0.999).mean())))
pd.DataFrame(qc).to_csv(COV / "pacbio_read_qc.tsv", sep="\t", index=False)
print("done")
