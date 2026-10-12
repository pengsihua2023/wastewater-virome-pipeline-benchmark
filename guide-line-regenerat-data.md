# Regenerating the simulated data: operating guide

**Current status (2026-09-03): all six datasets, Illumina and PacBio, have been generated and have passed acceptance.** The measured results are in the tables of the two gates in step 7. The raw files are not kept in version control (1.8–5.7 GB per file) and are listed in `.gitignore`; traceable intermediate products such as per-member counts and plan tables are kept under `regen/`.

**Goal:** generate six datasets (3× Illumina + 3× PacBio HiFi) in which **viral reads make up 4.0% of all reads** (phages 8.0%, bacteria and archaea 88.0%), the within-group relative abundances of the 708 members keep the designed proportions, no member is raised by a depth floor, and no member is silently lost.

**Method:** measurement and correction rather than derivation. How the simulators convert abundance into reads cannot be determined from the existing files, but the conversion is **linear** — the number of reads of a genome is proportional to the abundance given to it. So one small-scale trial run to measure the actual output of every member, followed by one linear correction, is enough to hit the target. This guide is designed accordingly.

**Scope:** this guide is written for an operator with access to the genome FASTA files, CAMISIM and PBSIM3 on the HPC. Every step that requires measurement is annotated with the exact commands and acceptance criteria.

---

## Step 0: starting point and prerequisites

0.1 is settled and can be adopted directly; one item in 0.2 is still missing.

### 0.1 Target composition: 88 / 8 / 4 refers to read share (settled; no need to ask the authors again)

The two seemingly contradictory sets of figures in the design document are two representations of the same numbers:

| Source | Bacteria/archaea | Phages | Viruses |
|---|---|---|---|
| `Simulation design reads` column (measured read counts of the source sample HTP-2024-01-15) | 5.62% | 1.00% | 93.38% |
| `Relative abundance in simulation (%)` column | 88.0% | 8.0% | **4.0%** |

**Basis of the decision.** We tested how the "relative abundance" column was computed from the "design reads"; fits of the two hypotheses:

| Group | Hypothesis A: abundance ∝ design reads (read share) | Hypothesis B: abundance ∝ design reads / genome length (copy share) |
|---|---|---|
| Bacteria/archaea | **r = 1.000000**, max deviation 8×10⁻⁹ | r = 0.080 |
| Phages | **r = 1.000000**, max deviation 6×10⁻⁸ | r = 0.995 |
| Viruses | **r = 1.000000**, max deviation 1×10⁻⁷ | r = 0.847 |

The abundance column equals the design reads normalised within each block and then scaled across the three blocks to 88/8/4, an exact identity (deviations at the
floating-point level). No step **divides by genome length**.

**Conclusion: 88 / 8 / 4 are the target read shares.** 93.38% is the raw read share of the source sample and 88/8/4 is the target after rescaling; the two are not
contradictory. This guide proceeds on the read-share interpretation.

**This decision determines how each platform is fixed:**

| | Abundance semantics | Matches the design? | Action |
|---|---|---|---|
| CAMISIM (Illumina) | copy share | ✗ deviates | **Must be fixed**: copy share = target read share ÷ genome length (step 3) |
| `simulate_pacbio.py` old version (2026-06-01) | base share | ✓ close to design | **Keep** the depth formula; fix read length and error handling (step 2) |
| `simulate_pacbio.py` corrected version (2026-06-02) | copy share | ✗ deviates | **Do not adopt its depth formula** — it pulls PacBio to the wrong side, the same as CAMISIM |

See `simulation-code-versions/README.md` for details.

### 0.2 Retrieve the missing materials

All of the following are required; none is currently in the data repository:

1. **CAMISIM configuration file** (`.ini`) — to confirm which abundance file it reads and whether it normalises internally. This is currently the only key item still missing
2. **Complete logs of both simulations** — to confirm the scale of the failures in step 2

The genome FASTA files are on the local machine (`CAMISIM-Generate/genomes/`, all 708 available and checked entry by entry against `genome_locations.tsv`),
and the four historical versions of `simulate_pacbio.py` are archived in `simulation-code-versions/`.

---

## Step 1: measure the genome lengths (mandatory, cannot be skipped)

The `est_genome_bp` column in the existing audit table `10_simulation_config_audit.xlsx` **is an estimate and some values are clearly wrong**: *Kobuvirus aichi* is listed as 752 bp (the real Aichivirus is about 8.2 kb) and *Sapporo virus* as 1,021 bp (real: about 7.5 kb). For 101 of the 151 viruses the "within-block median" was used directly.

**Any correction value computed from these lengths is wrong.** Lengths must be measured from the FASTA files:

```bash
python3 - <<'EOF'
import pandas as pd
from Bio import SeqIO

loc = pd.read_csv("genome_locations.tsv", sep="\t", header=None, names=["gid", "path"])
rows = []
for gid, path in zip(loc.gid, loc.path):
    try:
        L = sum(len(r.seq) for r in SeqIO.parse(path, "fasta"))
    except Exception as e:
        L = 0
        print(f"read failed {gid}: {path} -- {e}")
    rows.append({"gid": gid, "path": path, "length_bp": L})
g = pd.DataFrame(rows)
g["dom"] = g.path.str.extract(r"genomes/([^/]+)/")[0]
bad = g[g.length_bp == 0]
if len(bad):
    raise SystemExit(f"{len(bad)} genomes could not be read; resolve this before continuing:\n{bad.gid.tolist()}")
g.to_csv("genome_lengths.tsv", sep="\t", index=False)
print(g.groupby("dom").length_bp.describe()[["count","min","50%","max"]].to_string())
EOF
```

**Acceptance:** 708 rows and no member with length 0. Keep `genome_lengths.tsv` permanently and add this column to the design document.

---

## Step 2: the real PacBio defect is read-name collisions (determined by measurement)

### Read length is not the problem — this conclusion has been refuted experimentally

An earlier version of this guide stated that `--mean_length 15000` is incompatible with small genomes, on the basis that the audit table's
`obs/pred` ratio falls as genomes get smaller (phages 0.009, viruses 0.292). **That judgement was wrong.**

Measured diagnostic (6 genomes × 5 read lengths, 20× depth, counted as unique ZMWs in the BAM):

| Genome | Length | 15000 | 8000 | 4000 | 2000 | 1000 |
|---|---|---|---|---|---|---|
| genome_00083 | 910 bp | **1.000** | 1.000 | 1.000 | 1.000 | 1.050 |
| genome_00049 | 3,925 bp | **1.000** | 1.000 | 1.050 | 1.000 | 1.013 |
| genome_00002 | 7,600 bp | **1.000** | 1.050 | 0.947 | 1.000 | 0.987 |
| genome_00703 | 15,222 bp | 1.100 | 1.026 | 1.013 | 0.987 | 0.987 |
| genome_00174 | 43 kb | 1.000 | 1.000 | 0.968 | 0.972 | 0.982 |
| genome_00140 | 3.2 Mb | 0.975 | 0.973 | 0.976 | 0.976 | — |

The ratio is 0.95–1.10 for all sizes and all read lengths. **PBSIM3 truncates reads to the template length rather than discarding them**,
so a 910 bp genome still yields the full amount at a 15 kb read length. The read-length parameter does not need to be changed.

The size-dependent ratio in the audit table is a **consequence** of the defect described below, not its cause.

### The real defect: read-name collisions across genomes

The old script ran PBSIM3 separately for each genome, and PBSIM3 ZMW numbering **restarts at 1 every time**.
Measured evidence from the old data `data/Pacbio/pacbio_hifi_reads_sample_0.fastq.gz`:

- read names look like `@S1/1/ccs`; the largest ZMW number is only **4,339**, although the file holds 485,917 reads
- **8,108** duplicated names among the first 400,000 reads

When files were merged, reads with the same name overwrote each other, and the ones sacrificed were exactly those with small numbers — **the reads produced by small genomes**.
The consequences agree exactly with the measurements:

| | Expected under the old script formula | Actual in the old data |
|---|---|---|
| Bacteria/archaea | 68.5% | 99.8% |
| Phages | 24.1% | 0.09% |
| Viruses | 7.4% | 0.12% |

This also explains the anomalies in the assembly steps of nf-core/mag and KrakenMetaReads-nf on PacBio: **duplicate read names
cause Flye to abort**, which is exactly the problem pointed out in the comments of the 2026-06-02 corrected script. Attributing those anomalies
to "input contigs" was incomplete and they must be re-evaluated after regeneration.

### Two changes that must be made

**First, globally unique read names** (adopt the corrected version's approach, or add a genome prefix):

```bash
samtools fastq "${g}.ccs.bam" \
  | awk -v g="$gid" 'NR%4==1{sub(/^@/,"@"g"_"); print; next}{print}' > "${gid}.fq"
# yields @genome_00001_S1/1/ccs, which no longer collides across genomes
```

**Second, failures must no longer be skipped silently**: when PBSIM3 returned non-zero the old version only logged a warning and did `continue`,
so a genome was lost while the job still exited successfully. Record a list of failures and exit with a non-zero code.

### The depth floor must still be removed

`--min_depth 1.0 --min_reads 5` raised 64.9% of the viral members and compressed a 1.4×10⁵-fold span of designed abundance
into an observed 8-fold span. **Change it to `--min_depth 0 --min_reads 0`**; otherwise the linear relationship is broken and
the correction in step 5 fails.

### Performance note

With `--pass-num > 1`, PBSIM3 outputs **a multi-pass BAM rather than FASTQ**, and `ccs` is needed to produce the HiFi consensus.
If ccs is called once per contig, process start-up overhead dominates the run time (measured with 19,135 calls:
only 4 genomes advanced in 7 minutes). **Merge the several BAMs of one genome with `samtools cat` first, and call ccs once.**

---

## Step 3: generate the initial abundance files

The two platforms use different abundance semantics and **must be handled separately**:

- **CAMISIM (Illumina)** reads abundance as the **genome copy share**. Confirmed: back-calculating the current 88/8/4 design with the measured lengths predicts 0.0198% viral reads, measured 0.0200%, in agreement.
- **`simulate_pacbio.py` (PacBio)** reads abundance as the **share of total bases** (`allocated_bases = total_bases × abundance`).

```python
import pandas as pd

TARGET = {"Bacteria_Archaea": 88.0, "Bacteriophage": 8.0, "Viruses": 4.0}   # confirmed in step 0.1
g = pd.read_csv("genome_lengths.tsv", sep="\t")          # measured lengths from step 1

for i in range(3):
    d = pd.read_csv(f"distribution_{i}.txt", sep="\t", header=None, names=["gid", "ab"])
    m = g.merge(d, on="gid", how="left")
    assert m.ab.notna().all(), "some members are missing from distribution"

    # keep the original design proportions within a group; reallocate between groups by target read share
    share = m.groupby("dom").ab.transform(lambda s: s / s.sum()) * m.dom.map(TARGET)

    # Illumina: copy share ∝ target read share ÷ genome length
    ill = share / m.length_bp
    ill = ill / ill.sum()
    m.assign(ab=ill)[["gid", "ab"]].to_csv(
        f"distribution_{i}_illumina.txt", sep="\t", header=False, index=False,
        float_format="%.12f")

    # PacBio: the old script allocates by base share (bases ∝ abundance), and reads = bases / effective read length,
    #         so reads ∝ abundance / effective read length. For reads ∝ target share,
    #         abundance must be ∝ target share × effective read length. Step 5 corrects this further with measured results.
    #         Note: this formula presupposes the old depth formula and does not apply if the corrected version's formula is used.
    eff = m.length_bp.clip(upper=15000)      # change here if a different read-length strategy is chosen after step 2
    pb = share * eff
    pb = pb / pb.sum()
    m.assign(ab=pb)[["gid", "ab"]].to_csv(
        f"distribution_{i}_pacbio.txt", sep="\t", header=False, index=False,
        float_format="%.12f")

    # self-check: back-solve the Illumina read share
    chk = ill * m.length_bp
    print(f"sample {i} expected Illumina read share:",
          (100 * chk.groupby(m.dom).sum() / chk.sum()).round(2).to_dict())
```

**Acceptance:** the self-check output should be `{'Bacteria_Archaea': 88.0, 'Bacteriophage': 8.0, 'Viruses': 4.0}`.

> **Do not use `analysis_708/interim/abundance_correction_required.tsv`.** The correction values in that table = target ÷ measured output ratio, and the measured PacBio output ratio is itself depressed by the step-2 defect. If that table were adopted at the same time as the read-length problem is fixed, phages would be over-estimated by about 113-fold and viruses by about 3.4-fold.

---

## Step 4: first pass (full scale)

Abundance and output are linearly related, so **two passes are enough to hit the target**: the first pass uses the initial values of step 3 and measures the actual output of each member; the second pass uses the corrected abundances and gives the final data.

**The first pass must be run at full scale; do not scale it down.** A numerical test (8 independent repeats) shows:

| First-pass scale | Viral share after correction | Largest deviation over the three groups |
|---|---|---|
| 1/10 | 3.83 ± 0.14% | 1.57 percentage points |
| **Full scale** | **4.01 ± 0.02%** | **0.14 percentage points** |

Scaling down puts the read counts of low-abundance members in the Poisson-noise regime and distorts the correction factors.

```bash
# Illumina: CAMISIM config points to distribution_0_illumina.txt, with the production scale
camisim.py -c config_pass1.ini

# PacBio
python simulate_pacbio.py \
  --genome_locations genome_locations.tsv \
  --distribution     distribution_0_pacbio.txt \
  --output_dir       pass1_pacbio \
  --total_bases      5e9 \
  --method qshmm --model "$MODEL" --pass_num 10 \
  --mean_length <value determined in step 2> --sd_length <mean/7> --min_len 500 \
  --min_depth 0 --min_reads 0 \
  --seed 42 --threads 32
```

**Note `--min_depth 0 --min_reads 0`.** The original setting `--min_reads 5` raised 64.9% of the viral members and compressed a 1.4×10⁵-fold span of designed abundance into an observed 8-fold span; the 13/13 sentinel-virus detection on PacBio came mainly from this floor rather than from sensitivity. Keeping any floor breaks the linear relationship and invalidates the step-5 correction.

---

## Step 5: measure and correct

### 5.1 Per-member counts

```bash
# PacBio: each genome is simulated separately, so counts are direct
for f in pass1_pacbio/tmp/*.fastq.gz; do
  gid=$(basename "$f" | sed 's/_.*//')
  echo -e "${gid}\t$(zcat "$f" | awk 'END{print NR/4}')"
done > pass1_reads_pacbio.tsv

# Illumina: the CAMISIM reads mapping file can give per-genome origin;
#           if not available, fall back to group-level calibration (see 5.3)
```

### 5.2 Direct ratio correction

For each member: **new abundance = old abundance × target reads ÷ measured reads**. This is an exact linear back-solve and does not depend on any assumption about the simulator's internal behaviour.

```python
import numpy as np, pandas as pd

R_TOTAL = 16_642_947                    # actual total reads of the first pass (fill in the measured value)
TARGET  = {"Bacteria_Archaea": 88.0, "Bacteriophage": 8.0, "Viruses": 4.0}
MIN_CNT = 20                            # members below this count use the median factor of their size bin

g   = pd.read_csv("genome_lengths.tsv", sep="\t")
d   = pd.read_csv("distribution_0_pacbio.txt", sep="\t", header=None, names=["gid","ab"])
obs = pd.read_csv("pass1_reads_pacbio.tsv",    sep="\t", header=None, names=["gid","reads"])
m   = g.merge(d, on="gid").merge(obs, on="gid", how="left").fillna({"reads": 0})

share  = m.groupby("dom").ab.transform(lambda s: s / s.sum()) * m.dom.map(TARGET) / 100
m["target_reads"] = R_TOTAL * share

# only "sufficient target but nothing produced" counts as a simulation failure; go back to step 2 rather than hiding it by raising abundance
failed = m[(m.target_reads >= 10) & (m.reads == 0)]
if len(failed):
    raise SystemExit(
        f"{len(failed)} members have target ≥10 but output 0 -- simulation failure; go back to step 2 and adjust the read-length strategy:\n"
        f"{failed.gid.tolist()[:20]}")

# members with sufficient counts use their own ratio; the others use the median ratio of "same group × same size bin"
ratio = np.where(m.reads >= MIN_CNT, m.target_reads / m.reads.clip(lower=1), np.nan)
r = pd.Series(ratio)
bins = pd.cut(m.length_bp, [0, 2e3, 5e3, 1e4, 1.5e4, 2e4, 5e4, 1e5, 1e6, 1e10])
r = r.fillna(r.groupby([m.dom, bins], observed=True).transform("median")).fillna(r.median())

ab_new = m.ab * r
ab_new = ab_new / ab_new.sum()
m.assign(ab=ab_new)[["gid","ab"]].to_csv(
    "distribution_0_pacbio_final.txt", sep="\t", header=False,
    index=False, float_format="%.12f")

print(f"members with fewer than {MIN_CNT} reads that use the bin median factor: {(m.reads < MIN_CNT).sum()}/708")
```

Repeat steps 4–5 for samples 1 and 2 (each with its own distribution file and random seed).

### 5.3 If per-genome counts are not available for Illumina

Fall back to group-level calibration using the group shares obtained from Kraken2:

```python
k = m.groupby("dom").apply(lambda x: x.reads.sum() / x.target_reads.sum())
ab_new = m.ab / m.dom.map(k)      # scale proportionally within a group; within-group proportions are unchanged
ab_new = ab_new / ab_new.sum()
```

Group-level calibration is enough to guarantee the 88/8/4 target itself (it is a group-level metric), but the within-group distribution still depends on the accuracy of the analytical correction in step 3. **Prefer obtaining per-genome counts.**

---

## Step 6: second pass (final data)

Use `distribution_{0,1,2}_*_final.txt` at the same scale as the first pass, with different random seeds for the three samples (e.g. 42, 43, 44). The output of this pass is the final dataset.

**The total sequencing amount does not need to be increased.** The current depth is already 2.8 times what the rarest virus (designed abundance 1.7×10⁻⁵%) needs to reach an expected ≥1 read — the bottleneck was always the configuration, not the depth.

---

## Step 7: acceptance (all three gates must pass before analysis)

**Gate 1 is decisive and must be done first.** It parses the FASTQ headers directly and counts per member without any classifier.
`Verification-simulated-data-vs-design.md` showed that this method is feasible and complete: all 16,642,909 reads of the existing data
map to the 708 members, with 0 unmapped. The earlier Kraken2-dependent check
underestimated the severity of the defect (Kraken2 reported 3,326 viral reads where the raw data had 1,509).

### Measured results (2026-09-03, all six datasets passed both gates)

| Illumina | Bacteria/archaea | Phages | Viruses | Viruses with zero output |
|---|---|---|---|---|
| sample_0 | 88.138% | 7.876% | 3.986% | 0/151 |
| sample_1 | 88.063% | 7.920% | 4.017% | 0/151 |
| sample_2 | 88.254% | 7.875% | 3.870% | 0/151 |
| Target | 88.0% | 8.0% | 4.0% | — |

| PacBio HiFi | Bacteria/archaea | Phages | Viruses | Viruses with zero output | Whole community with zero output | Target ZMW achievement |
|---|---|---|---|---|---|---|
| sample_0 | 88.184% | 7.872% | 3.944% | 0/151 | 0/708 | 98.5% |
| sample_1 | 88.193% | 7.859% | 3.948% | 0/151 | 0/708 | 98.7% |
| sample_2 | 88.293% | 7.871% | 3.836% | 0/151 | 0/708 | 98.5% |
| Target | 88.0% | 8.0% | 4.0% | — | — | — |

Maximum deviation: Illumina 0.254 pp, PacBio 0.293 pp. The raw FASTQ files are not kept in version control (`.gitignore`
lists `regen/illumina/` and `regen/pacbio/*.fastq.gz`); per-member counts are retained in
`regen/per_member_s{0,1,2}.tsv` (Illumina) and `regen/pacbio/counts_s{0,1,2}.tsv`
(PacBio) for traceable re-checking.

On the PacBio side, two further independent defects were located and fixed along the way (see
`Verification-simulated-data-vs-design.md` §7 and
`simulation-code-versions/README.md`):

1. **Wrong abundance plan for multi-contig genomes**: `eff` had used the total genome length instead of the length of a single contig,
   which does not match PBSIM3 counting each contig independently (members with 471 contigs were over-produced 5-fold).
   It was changed to sum `max(1, contig length/15000)` over contigs as the coverage unit
   (`regen/cov_units.tsv`).
2. **Read-name collapse within multi-contig genomes**: PBSIM3 numbers each contig independently, the movie name
   is fixed and the hole number restarts at 1 for every contig, so merging directly makes reads of different contigs share a
   QNAME, and `samtools fastq` keeps only the first of reads with the same name (one measured case with 471 contigs collapsed
   from 29,409 to 64 reads). Before merging, the hole numbers of the contigs are now accumulated with offsets so that names are globally unique.

Eight members (across the 3 samples) with output more than 30% below target were caught by the automatic check: 6 are edge cases not covered by these two fixes (2 are genomes whose
ccs pass rate is inherently low under the error model, e.g.
*genome_00208* at only 26.8%), and 2 are ultra-low-abundance viruses with target=1 (ccs rejects them at random). All were re-run and repaired by
raising the depth multiplier (computed from each member's own measured pass rate, not a blanket factor); after the repair all reached their target counts.

### Gate 1: member-by-member check of FASTQ against design values (Illumina)

```bash
# 1a. count per source accession (CAMISIM encodes the source in the read header, e.g. @FJ755405.1-1735/1)
zcat <new data>/sample_*_01.fq.gz \
  | awk 'NR%4==1{h=substr($1,2); n=index(h,"-"); print (n>0 ? substr(h,1,n-1) : h)}' \
  | sort | uniq -c | sort -rn > acc_counts.txt
```

```python
# 1b. map to members and compare with the design values
import pandas as pd, gzip, os
BASE = "<root directory of the genome FASTA files>"
loc = pd.read_csv("genome_locations.tsv", sep="\t", header=None, names=["gid","path"])
loc["rel"] = loc.path.str.replace(r"^.*/genomes/", "", regex=True)
loc["dom"] = loc.rel.str.split("/").str[0]
loc["taxid"] = loc.rel.str.split("/").str[-1].str.split("_").str[0]

rows = []
for gid, rel, dom, tx in zip(loc.gid, loc.rel, loc.dom, loc.taxid):
    fp = os.path.join(BASE, rel)
    op = gzip.open if fp.endswith(".gz") else open
    with op(fp, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                rows.append((line[1:].split()[0], gid, dom, tx))
m = pd.DataFrame(rows, columns=["acc","gid","dom","taxid"])

cnt = [(p[1], int(p[0])) for p in (l.split() for l in open("acc_counts.txt")) if len(p) == 2]
c = pd.DataFrame(cnt, columns=["acc","reads"])
j = m.merge(c, on="acc", how="left").fillna({"reads": 0})

# completeness: must be 100% mapped
unmapped = c.reads.sum() - j.reads.sum()
assert unmapped == 0, f"{unmapped:,.0f} reads could not be mapped to members; check that the FASTA files and genome_locations are consistent"

per = j.groupby(["gid","dom","taxid"], as_index=False).reads.sum()
dom = per.groupby("dom").reads.sum()
print("group read share %:", (100*dom/dom.sum()).round(4).to_dict())
v = per[per.dom == "Viruses"]
print(f"viruses: 0 reads {int((v.reads==0).sum())}/151 | ≥1 read {int((v.reads>=1).sum())}/151 | "
      f"≥10 reads {int((v.reads>=10).sum())}/151 | median {v.reads.median():.0f}")
per.to_csv("fastq_reads_per_member.tsv", sep="\t", index=False)
```

| Metric | Pass criterion | Existing (defective) data |
|---|---|---|
| Read mapping completeness | **100%, 0 unmapped** | 100% (this item already passed) |
| Viral read share | **4.0% ± 0.4%** | **0.009%** ✗ |
| Phage read share | 8.0% ± 0.8% | 0.010% ✗ |
| Bacteria/archaea read share | 88.0% ± 2% | 99.98% ✗ |
| **Number of viruses with zero output** | **0 / 151** | **97 / 151** ✗ |
| Median viral read count | ≥ 50 reads | **0** ✗ |
| Per-member measured/design (scaled) | all within **0.5–2.0** | all about 0.0001 ✗ |

**Failing any one of these means the data must not enter analysis.** In particular the "number of viruses with zero output" must be exactly 0 — if any virus is completely absent, recall with 151 as the denominator is wrong.

### Gate 2: PacBio member-by-member check

The merged PBSIM3 read headers (`@S1/1/ccs`) contain no source information. The method actually adopted is to rewrite the read name to `@<gid>_S1/<hole>/ccs` already at the
`samtools fastq` output (`gid` being the internal identifier
`genome_00001`–`genome_00708`), so that the final merged FASTQ itself can be traced back to its source read by read, with no need to keep per-genome intermediate files:

```bash
samtools fastq "$w/ccs.bam" \
  | awk -v g="$gid" 'NR%4==1{sub(/^@/,"@"g"_");print;next}{print}' > "${gid}.fq"

# acceptance: count the per-genome fq files directly (generation already writes one file per genome, so no need to split the merged file again)
for f in genome_*.fq; do
  echo -e "$(basename "$f" .fq)\t$(( $(wc -l < "$f") / 4 ))"
done > pacbio_reads_per_member.tsv
```

The pass criteria are the same as for gate 1. The measured results for the three samples are in the "Measured results" table above.

### Gate 3: script-level re-check

```bash
python3 analysis_708/scripts/05_realised_composition.py      # composition as seen by the classifier
python3 analysis_708/scripts/08_simulation_config_audit.py   # depth floor and output ratio
python3 analysis_708/scripts/07_designed_vs_observed_reads.py
```

| Metric | Pass criterion | Existing data |
|---|---|---|
| Number of members with `floored` True | **0** | viruses 98/151 |
| Median `pb_ratio_obs_over_pred` per group | 0.9–1.1 | bacteria 1.475, viruses 0.292, phages 0.009 |

Gate 3 is **auxiliary**: it goes through a classifier and underestimates the defect (Kraken2 misclassifies bacterial reads as viral and inflates the viral read count). **Gates 1 and 2 take precedence.**

### After all three gates pass

```bash
for s in analysis_708/scripts/[0-9][0-9]_*.py; do python3 "$s" || break; done
```

`13_consistency_check.py` must report **0 failures**.

---

## Step 8: what the manuscript must re-measure after regeneration

| Conclusion | Re-measurement needed? | Notes |
|---|---|---|
| Absolute recall / precision / F1 of each pipeline | **Must be re-measured** | Of the 151 denominators, 97 viruses are not in the data at all (0 reads). With the 54 actually present as the denominator, the true recall of the leading pipelines is 0.80–0.83 rather than 0.30 — the current values seriously **underestimate** pipeline capability |
| Ranking of the seven pipelines | **Must be re-measured** | The ranking changes as soon as the 54 detectable viruses are used as the denominator: sourmash rises to first (0.833), GOTTCHA2 falls from first to fifth (0.667) |
| "Detected from nothing" counted as true positives | **Must be removed** | Some viruses counted as TP have 0 reads in the data and cannot really have been detected (CLARK 7, taxprofiler 2). The current scoring inflates both recall and precision |
| Sentinel detection rates (3/13, 13/13) | **Must be re-measured** | Part of the 3/13 in Illumina is due to the sentinel viruses being absent from the data; the 13/13 in PacBio comes mainly from the depth floor |
| Threshold filtering lowers F1 | **Will necessarily change** | With tens to hundreds of reads per virus, the ≥10 threshold no longer deletes the signal |
| reads vs contigs control (assembly loses viral signal) | Expected to hold | A within-data control independent of absolute depth, but the magnitude will change |
| ±fastp control (false positives double) | Expected to hold | Same as above |
| Three pipelines return zero on PacBio | Expected to hold | A qualitative fact; but GOTTCHA2 additionally has an unvalidated patch problem |

The extent of the defect in the existing data is in `Verification-simulated-data-vs-design.md`: of 151 viruses, 97
have not one read in the raw FASTQ, and the total viral share is 0.009%. The absolute metrics currently in the manuscript are therefore neither a
"lower bound" nor a "conservative estimate" but **numbers measured on a community that does not match the design**, and must not be used to
characterise pipeline capability.

After regeneration, the viral reads per Illumina library will rise from 3,326 to about 666,000, the number of viruses with an expected <1 read falls from **94/151 to 0/151**, and the median expected reads per virus rises from 0.47 to 93.

---

## Appendix: problems unrelated to the simulation that also block submission

Regenerating the data does not solve the following:

1. On Illumina the Kraken2 reports of nf-core/taxtriage and nf-core/taxprofiler are **byte-identical** (MD5-verified) — confirm whether taxtriage was really run independently
2. On PacBio the three read-level reports of nf-core/mag are **byte-identical to each other** — they must be re-run per sample
3. The GOTTCHA2 PacBio patch has **not been validated with a positive control**, so its zero recall cannot be distinguished from a patch defect
4. Version information is missing: CAMISIM and the read simulator it calls, GOTTCHA2 and CLARK versions and database build dates, the exact taxtriage revision, Nextflow and container runtime versions
5. `databases.csv` labels the taxprofiler database as `Viral_ref`, but the actual path points to the Kraken2 **Standard** library

---

## Confidence in this guide

| Content | Basis | Confidence |
|---|---|---|
| CAMISIM interprets abundance as copy number | Back-calculation with measured lengths: 0.0198% vs measured 0.0200% | **Confirmed** |
| Copy share = target read share ÷ genome length | Back-solve check returns to 88/8/4 | **Confirmed** |
| The depth floor raises 64.9% of viral members | Recomputed from the script formula | **Confirmed** |
| Small genomes are strongly under-produced, monotonically with size | Measured/predicted ratio of the 708 members | **Confirmed** |
| The exact mechanism of the under-production (discard vs failure) | Not tested, logs missing | **To be determined by measurement in step 2** |
| The read-length parameter values after step 2 | Determined by experiment; this guide does not presuppose them | **To be determined by the operator** |
| `est_genome_bp` in the audit table | Estimate, some values clearly wrong | **Unusable; must be re-measured** |

The guide adopts "measure and correct" rather than "analytical solution" precisely so that the last two uncertainties do not affect the result: however the simulator converts internally, steps 4–5 measure it and correct it away.
