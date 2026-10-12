# Versions of the PacBio simulation script

This directory collects all versions of `simulate_pacbio.py` found on the local machine, together with
two versions of the run script, in order to establish **which version generated the data used in the analysis**.

## Version list

| File | Original location | Date | Size (bytes) | Notes |
|---|---|---|---|---|
| `simulate_pacbio_v2026-03-18_earliest.py` | `CAMISIM/` | 2026-03-18 | 10,763 | Earliest version |
| `simulate_pacbio_v2026-05-29.py` | `CAMISIM-Generate/` | 2026-05-29 | 14,504 | |
| `simulate_pacbio_v2026-06-02_FIXED.py` | `CAMISIM-Generate/generation-code/` | 2026-06-02 | 17,412 | Corrected version, **not used to generate the existing data** |
| `run_pacbio_sample_v2026-06-07.sh` | `CAMISIM-Generate/generation-code/` | 2026-06-07 | 1,235 | Companion of the corrected version |

The 2026-06-01 script (and its run script) generated the first PacBio dataset, which was found to be defective (see below) and was discarded; these two files have been removed from the repository. The data used in the paper were regenerated (`data/regen/`, `docs/guide-line-regenerat-data.md`).

## How we established that the data came from the old version

The two versions write different read header formats:

| Version | Code | Output header |
|---|---|---|
| Old | `out.write(header)` | `@S1/1/ccs np=10 rq=0.998301` |
| Corrected | `out.write(f"@r{total_reads} {orig_id}\n")` | `@r1 S1/1/ccs np=10 rq=0.998301` |

The headers in `data/Pacbio/pacbio_hifi_reads_sample_0.fastq.gz` look like `@S1/1/ccs …`,
**without the `@r{N}` prefix, so the file was produced by the old version**.

## What the corrected version changed

### 1. Depth allocation formula

Original comment in the corrected version:

> the previous version used `depth = total_bases * abundance / size`, which makes
> reads ∝ abundance (ignoring genome size). That massively over-represented
> small-genome organisms (e.g. viruses) relative to CAMISIM — it inflated the
> viral read fraction by ~600x.

| | Formula | Effect |
|---|---|---|
| Old | `depth = total_bases × abundance / genome_size` | bases ∝ abundance, **independent of genome length** |
| Corrected | `depth = total_bases × abundance / Σ(abundance × size)` | bases ∝ abundance × genome length, **consistent with CAMISIM** |

### 2. Unique read IDs — the actual root cause of the defect in the old PacBio data

PBSIM3 is run separately for each genome, and ZMW numbering restarts at 1 every time, so IDs are duplicated after merging.
The corrected version assigns every read a globally unique ID and keeps the original ID as a comment.

**Measured evidence** (`data/Pacbio/pacbio_hifi_reads_sample_0.fastq.gz`): among 485,917 reads
the largest ZMW number is only **4,339**, and the first 400,000 reads contain **8,108** duplicated names. When files were merged,
reads with the same name overwrote each other, and the losers were the genomes with small ZMW numbers — **small genomes that
produced few reads**. Consequence:

| Group | Expected under the old formula | Actual in the old data |
|---|---|---|
| Bacteria/archaea | 68.5% | **99.8%** |
| Phages | 24.1% | **0.09%** |
| Viruses | 7.4% | **0.12%** |

In other words, the old abundance formula (above) was itself close to the design; it was the read-name collision that erased the viruses.
Duplicate names also **caused Flye to abort**, consistent with the anomalies in the assembly steps of nf-core/mag and KrakenMetaReads-nf on PacBio
— attributing those anomalies to "input contigs" alone was therefore incomplete.

## Verdict: the old version's abundance semantics were correct; the corrected version went in the wrong direction

The corrected version changes the PacBio depth formula to match CAMISIM. But **it is the CAMISIM side that deviates from the design**.

### Evidence

We tested how the "relative abundance" column of the design document was computed from the "design reads":

| Group | Hypothesis A: abundance ∝ design reads (read share) | Hypothesis B: abundance ∝ design reads / genome length (copy share) |
|---|---|---|
| Bacteria/archaea | **r = 1.000000**, max deviation 8×10⁻⁹ | r = 0.080 |
| Phages | **r = 1.000000**, max deviation 6×10⁻⁸ | r = 0.995 |
| Viruses | **r = 1.000000**, max deviation 1×10⁻⁷ | r = 0.847 |

The abundance column equals the design reads normalised within each block and then scaled to 88/8/4 across the three blocks — an exact identity.
No step divides by genome length. **Therefore 88/8/4 are target read shares, not copy shares.**

### Conclusion

| | Abundance semantics | Matches the design? |
|---|---|---|
| Old PacBio (bases ∝ abundance) | base/read share | **✓ close to design** |
| Corrected PacBio (bases ∝ abundance × length) | copy share | ✗ deviates |
| CAMISIM (Illumina) | copy share | ✗ deviates |

The comment in the corrected version says the old version "inflated the viral read fraction by ~600x" — **that factor of 600 is not an inflation: the old
version was close to the design and CAMISIM was too low**. Redoing the PacBio data with the corrected version would push the PacBio viral share from the current 0.18% down
to the 0.01% range of CAMISIM, which is worse.

**What really needs fixing is the CAMISIM input on the Illumina side**: set the copy abundance to
`target read share ÷ genome length` (see step 3 of `guide-line-regenerat-data.md`).

### Trade-offs when regenerating

| Source | What to adopt |
|---|---|
| **Old version (2026-06-01)** | Depth formula `depth = total_bases × abundance / genome_size` |
| **Corrected version (2026-06-02)** | Unique read IDs (`@r{N} {orig_id}`), so that Flye does not abort on duplicate IDs |
| Needed in both versions | PBSIM3 failures must no longer be skipped silently (the old version only logged a warning and then did `continue`, so a genome was lost while the job still exited successfully); remove the depth floor with `--min_depth 0 --min_reads 0` |

## Differences in run parameters

| Parameter | Old run script | 2026-06-07 version |
|---|---|---|
| `--total_bases` | 5e9 | **10e9** |
| `--min_depth` | 1.0 | **0** |
| `--min_reads` | 5 | **1** |
| `--mean_length` | 15000 | 15000 (unchanged) |

The depth floor was relaxed substantially in the new script (`--min_reads 5 → 1`), which shows that the problem had been noticed.

`--mean_length 15000` was not changed, and **it did not need to be**. It had been suspected that this read length is incompatible with genomes shorter than 15 kb
(132/151 viruses, 67/67 phages); **experiments ruled this out**: across 6 genomes (910 bp to 3.2 Mb)
× 5 read lengths (15000 to 1000) the observed/theoretical ZMW ratio was always 0.95–1.10, and 1.000 for the 910 bp genome at a 15 kb
read length. PBSIM3 truncates reads to the template length rather than discarding them. See
step 2 of `guide-line-regenerat-data.md`.
