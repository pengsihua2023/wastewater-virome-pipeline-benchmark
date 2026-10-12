# Verification report: member-by-member check of the simulated data against the design

**Object checked:** `data/Illumina/sample_708_0_01.fq.gz` (Illumina sample_0 used in this project's analysis, R1)
**Date of check:** 2026-08-18
**Conclusion: the simulated data were not generated according to the design. 97 of the 151 viruses are completely absent from the data (0 reads), and the total viral share is 0.009%, three to four orders of magnitude below the design.**

> **Update (2026-09-03): the data have been regenerated following the conclusions of this report, and all six datasets (3× Illumina + 3× PacBio)
> passed the member-by-member acceptance check.** Results and methods are in step 7 of `guide-line-regenerat-data.md`. This report documents the check of the
> **original defective data** and is kept as the historical basis for locating the problem; it does not describe the state of the data currently in use.

---

## 1. Why this check was needed

All earlier evaluations of the simulated data were based on the **Kraken2 classification report**, one step removed from the raw data. An alternative explanation therefore could not be excluded: *the FASTQ contains many viral reads, but Kraken2 failed to classify them correctly.* If that were true, the problem would lie in the classification step rather than the simulation step, and the remedy would be entirely different.

This check bypasses all classifiers and **parses the FASTQ read headers directly**. CAMISIM encodes the source sequence accession of every read in its header:

```
@FJ755405.1-1735/1
 └────┬────┘
 source accession
```

This allows an exact count of the reads produced by each genome, independent of any classification tool.

---

## 2. Methods

```bash
# count reads per accession
zcat data/Illumina/sample_708_0_01.fq.gz \
  | awk 'NR%4==1{h=substr($1,2); n=index(h,"-"); print (n>0 ? substr(h,1,n-1) : h)}' \
  | sort | uniq -c | sort -rn > acc_counts.txt
```

The accession → member mapping was built from the headers of the 708 genome FASTA files (local path `CAMISIM-Generate/genomes/`, matching `genome_locations.tsv` entry by entry): the 708 genomes contain 19,804 sequences in total.

**Mapping completeness: all 16,642,909 reads in the FASTQ map to the 708 members; 0 reads are unmapped.** The statistics below are therefore complete, with no omissions or ambiguous assignments.

---

## 3. Results: group level

| Group | Genomes | Actual reads | Actual share | Design share ᵃ | Abundance-column target ᵇ |
|---|---|---|---|---|---|
| Bacteria/archaea | 490 | 16,639,725 | **99.981%** | 5.62% | 88% |
| Phages | 67 | 1,675 | 0.010% | 1.00% | 8% |
| **Viruses** | **151** | **1,509** | **0.009%** | **93.38%** | **4%** |

ᵃ Obtained by normalising the `Simulation design reads` column of the design document (it reproduces the read-count distribution of the source sample HTP-2024-01-15)
ᵇ The `Relative abundance in simulation (%)` column of the design document

The two columns are not contradictory: the abundance column is the design reads normalised within each block and then scaled across the three blocks to 88/8/4, an exact identity
(fit r = 1.000000 for all three groups, maximum deviation 8×10⁻⁹). **88/8/4 are the target read shares**,
and 93.38% is the raw read share of the source sample. See step 0.1 of `guide-line-regenerat-data.md` for the derivation.

Measured against the designed target read share of 4%, the actual viral share is **440 times too low**; measured against the 93.38% of the source sample it is
about 10,000 times too low. Bacteria crowd out almost all of the reads.

---

## 4. Results: per virus

| Criterion | Viruses / 151 |
|---|---|
| No reads at all (0) | **97** |
| < 10 reads | 138 |
| ≥ 1 read | **54** |
| ≥ 2 reads | 33 |
| ≥ 5 reads | 20 |
| ≥ 10 reads | 13 |
| ≥ 100 reads | 3 |

The **median** viral read count is 0, and the maximum is 868 (*Mamastrovirus 1*).

Design values, scaled by library size, compared one by one:

| Virus | Design reads (scaled) | Actual in FASTQ | Achieved |
|---|---|---|---|
| *Mamastrovirus 1* | 9,150,348 | 868 | 0.01% |
| *Norwalk virus* | 1,360,004 | 143 | 0.01% |
| *Sapporo virus* | 1,154,437 | 119 | 0.01% |
| *Mamastrovirus 6* | 946,384 | 82 | 0.01% |
| *Mamastrovirus 9* | 345,109 | 31 | 0.01% |
| *Orthopicobirnavirus hominis* | 147,100 | **3** | 0.00% |

**The achieved fraction stops neatly at 0.01% for all viruses.** This uniformity is itself diagnostic: it is not a random failure of individual genomes but a single systematic unit-conversion error acting on all viruses at once.

The complete per-member data are in `analysis_708/interim/fastq_reads_per_member.tsv`.

---

## 5. Kraken2 is excluded, and the direction is opposite to the concern

| | Viral reads |
|---|---|
| Actually present in the raw FASTQ | **1,509** |
| Kraken2 (taxprofiler) report | 3,326 |

Kraken2 reports **2.2 times more** viral reads than are actually present. The excess consists of false positives from bacterial reads misclassified as viral. Therefore:

- **the classification step is not the source of the defect** and can be fully excluded
- the earlier Kraken2-based estimate of the viral share, 0.020%, was **optimistic**; the true value is 0.009%

---

## 6. Impact on the existing analysis conclusions

### 6.1 The denominator of recall is wrong

97 viruses have not a single read in the data, so no pipeline could possibly detect them. Computing recall over 151 counts "missing data" as "pipeline misses".

Taking Illumina sample_0 as an example, using the 54 viruses actually present in the data as the denominator:

| Pipeline | Reported in the manuscript (over 151) | Actually detectable (over 54) | Detected from nothing ᶜ |
|---|---|---|---|
| sourmash | 0.298 | **0.833** | 0 |
| nf-core/taxprofiler | 0.298 | **0.796** | 2 |
| nf-core/taxtriage | 0.298 | **0.796** | 2 |
| nf-core/mag | 0.298 | **0.796** | 2 |
| GOTTCHA2 | 0.245 | 0.667 | 1 |
| CLARK | 0.185 | 0.389 | **7** |
| KrakenMetaReads-nf | 0.126 | 0.352 | 0 |

ᶜ Number of viruses counted as true positives although they have 0 reads in the raw data

**Actual pipeline performance is far better than reported in the manuscript**: the true recall of the leading pipelines is about 0.80–0.83, not 0.30. The absolute metrics currently in the manuscript seriously underestimate pipeline capability.

### 6.2 Some "true positives" cannot be true

The viruses in the "detected from nothing" column have no reads in the data, yet the pipelines report them and, because they are in the ground-truth set, they are counted as true positives. These are in fact false positives; CLARK has 7 of them. The current scoring counts them as TP, inflating both recall and precision.

### 6.3 Rankings may change

Ranked by recall over 54, sourmash comes first (0.833), GOTTCHA2 drops to fifth (0.667) and CLARK to second-to-last (0.389), whereas in the manuscript GOTTCHA2 has the highest F1 on Illumina. The ranking conclusions must be re-measured on correct data.

---

## 7. PacBio: a different defect, but equally serious

**PacBio data cannot be checked member by member with the method of section 2.** The read headers written by PBSIM3 are `@S1/1/ccs`,
and the source information is discarded when files are merged. The headers themselves, however, expose another defect.

### 7.1 Read-name collisions across genomes

PBSIM3 ZMW numbering **restarts at 1 for every genome**, and the old script called
PBSIM3 once per genome for 708 genomes and merged the outputs directly without making names unique. Measured on `data/Pacbio/pacbio_hifi_reads_sample_0.fastq.gz`:

| Metric | Value |
|---|---|
| Total reads | 485,917 |
| Maximum ZMW number | **4,339** |
| Duplicated names among the first 400,000 reads | **8,108** |

When files were merged, reads with the same name overwrote each other, and the ones lost were those with small numbers — the **small genomes that produced few reads**.
The consequences agree with the measurements:

| Group | Expected under the old script formula | Actual in the old data |
|---|---|---|
| Bacteria/archaea | 68.5% | **99.8%** |
| Phages | 24.1% | **0.09%** |
| Viruses | 7.4% | **0.12%** |

Hence **the loss of viruses in Illumina and in PacBio has two different, unrelated causes**: in Illumina, CAMISIM misread the
read share as a copy share; in PacBio, read names collided. Each chain independently pushed viruses below 0.1%.

### 7.2 A hypothesis that has been refuted

It had been suspected that the defect came from `--mean_length 15000` being incompatible with small genomes (132/151 virus genomes are shorter than 15 kb).
**A direct experiment refuted this hypothesis**: across 6 genomes (910 bp to 3.2 Mb) × 5 read lengths (15000 to 1000),
the observed/theoretical ZMW ratio was always 0.95–1.10, and 1.000 for the 910 bp genome at a 15 kb read length.
PBSIM3 truncates reads to the template length rather than discarding them. The read-length parameter does not need to be changed.

### 7.3 Knock-on effects

Duplicate read names cause Flye to abort. The anomalies in the assembly steps of nf-core/mag and KrakenMetaReads-nf on PacBio
had been attributed to "input contigs"; that attribution is incomplete and must be re-evaluated after regeneration.

### 7.4 Scope of the check

**Only R1 of Illumina sample_0 was checked.** R2 and sample_1 and sample_2 were not checked, but the three samples were generated from
the same configuration and are expected to behave the same way.

---

## 8. How to reproduce

```bash
# 1. count reads per accession
zcat data/Illumina/sample_708_0_01.fq.gz \
  | awk 'NR%4==1{h=substr($1,2); n=index(h,"-"); print (n>0 ? substr(h,1,n-1) : h)}' \
  | sort | uniq -c | sort -rn > acc_counts.txt

# 2. build the accession -> member mapping (requires the genome FASTA files)
#    parse the '>' lines of the 708 FASTA files listed in genome_locations.tsv

# 3. merge counts and summarise by group and by member
#    result table: analysis_708/interim/fastq_reads_per_member.tsv
```

---

## 9. Conclusions and recommendations

1. **The simulated data have a major defect and cannot serve as the basis for a pipeline benchmark.** 97 of the 151 viruses are completely absent, and the viral share is three to four orders of magnitude below the design.
2. **The defect is in the simulation step; the classification step has been excluded.**
3. **The data must be regenerated.** The steps are in `guide-line-regenerat-data.md`.
4. **Acceptance after regeneration must use the method of this report** — parse the FASTQ headers member by member and compare with the design values, rather than relying on classifier output. Step 7 of the guide has been updated accordingly.
5. Until the new data are produced, none of the absolute metrics, pipeline rankings or sentinel detection rates in the manuscript should be released. Within-data controlled comparisons (reads vs contigs, ±fastp) do not depend on absolute depth and are expected to still hold.

**Outcome (2026-09-03):** all of the above recommendations have been implemented. The six datasets were regenerated and passed the member-by-member acceptance check,
with viral shares of 3.84–4.02% (three Illumina samples) and 3.84–3.95% (three PacBio samples); the number of the 151 viruses with zero output
is 0 in all six samples, consistent with the 88/8/4 target required by this report (maximum deviation 0.293 pp). Detailed data and the two further multi-contig-related defects found on the PacBio side
(wrong abundance plan, collapsed read names) are in `guide-line-regenerat-data.md`, step 7, "Measured results". **Recommendation 5 is now lifted**: the manuscript may re-measure absolute
metrics and pipeline rankings on the new data.
