# Data audit report — 708-member simulated wastewater benchmark

Prepared before any performance analysis was run. Every file under `Illumina/`
and `Pacbio/` was opened and classified; the complete list with sizes, MD5
checksums and a per-file decision is `01_file_inventory.xlsx` (302 files).
No input file was modified, moved or deleted.

---

## 1. What the project actually contains

| Item | Finding |
|---|---|
| Platforms | 2 — `Illumina/`, `Pacbio/` |
| Samples per platform | 3 (`sample_0`, `sample_1`, `sample_2`) — confirmed for **both** platforms |
| Pipelines | 7, present for both platforms: CLARK, GOTTCHA2, KrakenMetaReads-nf, mag-nf, sourmash, TaxProfiler, TaxTriage |
| Pipeline × platform × sample cells | 42; all 42 contain at least one parsable taxonomic table |
| Files inventoried | 302 |
| Files parsed into the detection table | 84 primary + sensitivity result tables |
| Files checked but not scored | 126 (redundant, derived or read-/contig-level files — reason given per file) |
| Context-only files | 95 (SLURM logs, submission scripts, Nextflow configs, sample sheets) |
| Raw FASTQ, databases, CAMISIM configuration | **absent from the project** |

### 1.1 The community is not 708 viruses

`all_species_list_for_simulation_data.docx` describes a **708-member microbial
community**, not 708 viruses. The composition extracted from its seven tables is:

| Section | n | Designed relative abundance |
|---|---|---|
| Bacteria / Archaea | 490 | 88.00% |
| Bacteriophage (ssRNA phage, all at NCBI rank S1) | 67 | 8.00% |
| Viruses — automated selection | 138 | 3.995% |
| Viruses — manually added (3 sub-categories) | 13 | 0.005% |
| **Total** | **708** | **100.00%** |

Verification: 708 rows, 708 unique names, 708 unique NCBI taxids, 708 unique
normalised match keys — **no duplicates, no synonyms, no unresolvable entries.**
Viral members total 218 (67 phage + 151 non-phage virus).

Any statement of the form "708 viruses" is incorrect and was not carried into
the analysis or the paper.

### 1.2 The document carries two mutually inconsistent specifications

`Simulation design reads` equals the actual Kraken2 clade reads of sample
HTP-2024-01-15 for every automatically selected member (verified: the two
columns are identical for all 695 automated entries). It is therefore a single,
community-wide **read-count** specification, and it implies the source sample's
own composition:

| Block | Designed reads | Share of the 708-taxon read budget |
|---|---|---|
| Bacteria/archaea | 57,147 | 5.62% |
| Bacteriophage | 10,214 | 1.00% |
| Virus | 949,995 | **93.38%** |

`Relative abundance in simulation (%)` is the column that was renormalised —
**within each block** — to 88 / 8 / 4. The ratio of design reads to relative
abundance is ~649 in the bacterial block, ~1,277 in the phage block and ~237,500
in the viral blocks, which is the signature of that per-block rescaling.

The two columns therefore specify different communities: 93.4% viral reads
versus 12% viral. Section 2 establishes which one reached the simulator.

---

## 2. The design's viral share was intended as reads but consumed as genome copies

Measured from the read-level Kraken2 report each platform's TaxProfiler run
produced against the same database (`05_realised_composition.py`):

| Platform | Sample | Total reads | % unclassified | % Bacteria | % Archaea | % Viruses |
|---|---|---|---|---|---|---|
| Illumina | sample_0 | 16,642,897 | 6.06 | 93.85 | 0.022 | **0.0200** |
| Illumina | sample_1 | 16,643,343 | 6.08 | 93.82 | 0.021 | **0.0201** |
| Illumina | sample_2 | 16,642,601 | 6.12 | 93.79 | 0.021 | **0.0198** |
| PacBio | sample_0 | 449,134 | 0.77 | 98.76 | 0.034 | **0.2220** |
| PacBio | sample_1 | 490,009 | 0.72 | 98.87 | 0.034 | **0.1829** |
| PacBio | sample_2 | 518,608 | 0.68 | 98.92 | 0.032 | **0.1849** |

The design document specifies **12%** of the community as viral (8% phage + 4%
virus), yet the measured viral **read** fraction is 0.020% (Illumina) and 0.20%
(PacBio). These two numbers are not the same quantity, and the gap between them
is explained — quantitatively — by genome length.

### 2.1 Direct comparison against the document's own read counts

The design document states a read count for every member, so the simplest test
is to compare it with what the datasets delivered
(`scripts/07_designed_vs_observed_reads.py`, `09_designed_vs_observed_reads.xlsx`).
Scaling the design read budget up to the realised Illumina library size
(16.6 M reads versus a 1.02 M design budget, a factor of 16.4):

| Member | Designed reads | Expected at realised library size | Observed (mean of 3) | Shortfall |
|---|---|---|---|---|
| *Mamastrovirus 1* | 559,347 | 9,150,369 | 606 | 15,100× |
| Norwalk virus | 83,135 | 1,360,007 | 135 | 10,100× |
| Sapporo virus | 70,569 | 1,154,440 | 12 | 93,600× |
| *Mamastrovirus 6* | 57,851 | 946,386 | 81 | 11,600× |
| *Kobuvirus aichi* | 13,393 | 219,096 | 2 | 131,000× |
| *Orthopicobirnavirus hominis* | 8,992 | 147,100 | 0 | not detected |

By block, reads delivered per designed read, normalised to bacteria/archaea:
phage 1/1,510 and virus 1/212,000 on Illumina; 1/201 and 1/14,900 on PacBio.

### 2.1.1 Which column reached the simulator

Weighting each candidate specification by genome length gives a decisive
discrimination against the measured viral read share:

| Hypothesised CAMISIM input | Predicted viral read share |
|---|---|
| `Simulation design reads` (viral 93.4% of reads) | 2.84% |
| `Relative abundance` 88 / 8 / 4 (viral 12%) | **0.017%** |
| **Measured, Illumina** | **0.017%** |

The measured value matches the relative-abundance hypothesis to three decimal
places and is 170-fold below the design-reads hypothesis. **CAMISIM was given the
88 / 8 / 4 relative-abundance column, and consumed it as genome copies.**

The viral read content therefore lost ground twice: the per-block renormalisation
cut it from 93.4% to 12% (a deliberate design choice, ~8×), and the copy-versus-
read interpretation cut it from 12% to 0.017% (~700×, unintended).

### 2.2 Diagnostic test

CAMISIM consumes an abundance file expressed as **genome copy number**, and the
number of reads a genome yields scales with copy number × genome length. If the
design percentages are copy abundances, then for each domain

    (total observed clade reads) / (total designed abundance)

normalised to bacteria should equal that domain's mean genome length divided by
the mean bacterial genome length — not 1. `scripts/06_abundance_scale_diagnosis.py`
measures this (mean of three samples per platform):

| Platform | Domain | Observed reads per unit abundance (bacteria = 1) | Observed deficit | Deficit predicted by genome length |
|---|---|---|---|---|
| Illumina | Bacteria/archaea | 1.000 | — | — |
| Illumina | Bacteriophage (~4 kb) | 0.0012 | **847×** | 1000× |
| Illumina | Virus (~7 kb) | 0.0016 | **645×** | 571× |
| PacBio | Bacteria/archaea | 1.000 | — | — |
| PacBio | Bacteriophage | 0.0091 | 110× | 1000× |
| PacBio | Virus | 0.0228 | 44× | 571× |

(Predicted values use approximate genome sizes — 4 Mb bacterial, 4 kb ssRNA
phage, 7 kb for a viral block dominated by small RNA viruses such as
*Mamastrovirus 1*, Norwalk virus and Sapporo virus — and are order-of-magnitude
expectations, not measurements.)

On Illumina the observed and predicted deficits agree to within ~20%. Feeding
the observed per-domain yields back through the design abundances predicts a
viral read share of 0.018%, against 0.021% measured. **The simulator therefore
executed its abundance file faithfully — but that file encoded genome copies,
not the read shares the design intended (Section 2.2).**

### 2.3 What the finding actually is

Both of the document's specifications are expressed in **reads**: the design-read
column literally states a read count per member, and the 88 / 8 / 4 split is a
read-share target (confirmed by the study authors). Neither was ever a
genome-copy specification. Note also that 88 / 8 / 4 is already a large,
deliberate downscaling of the source sample's own 5.6 / 1.0 / 93.4 read split,
not an inflation of it — the viral read target is conservative on both readings.

CAMISIM, however, consumes its abundance file as **genome copies**. Passing
88 / 8 / 4 to it therefore produced a viral *read* share of 0.018–0.020%
(Illumina) rather than 12%. The simulator behaved correctly for the input it was
given; the input encoded the wrong quantity.

**The datasets consequently miss the intended viral read composition by ~675×
on Illumina and ~65× on PacBio.** This is a data defect requiring re-simulation,
not merely a wording problem — although the design document's failure to
distinguish a copy share from a read share is what allowed it.

### 2.3.1 The correction required

Because read yield per unit copy abundance is proportional to genome length, the
copy abundances needed to deliver the intended read shares are
(`interim/abundance_correction_required.tsv`):

| Platform | Domain | Current design (%) | Required copy abundance (%) |
|---|---|---|---|
| Illumina | Bacteria/archaea | 88.0 | 0.93 |
| Illumina | Bacteriophage | 8.0 | 71.7 |
| Illumina | Virus | 4.0 | 27.3 |
| PacBio | Bacteria/archaea | 88.0 | 7.7 |
| PacBio | Bacteriophage | 8.0 | 77.0 |
| PacBio | Virus | 4.0 | 15.3 |

Two practical consequences. First, the required abundance file **differs between
platforms**, because read-length capping changes the yield per copy for small
genomes; a single abundance file cannot give both platforms a 12% viral read
share. Second, the general and platform-independent fix is to specify each
member's copy abundance as (target read share) ÷ (genome length), which requires
a genome-length column that the current design document does not carry.

### 2.4 Why the two platforms differ ~10-fold

The PacBio deficits (110× and 44×) are about 8–15× smaller than the Illumina
ones (847× and 645×). This is a read-length effect: a HiFi read approaching or
exceeding the length of a whole 7 kb viral genome means each viral molecule
yields roughly one read regardless of its length, whereas a 4 Mb bacterial
genome is still cut into hundreds of reads. Per unit copy abundance, viruses
therefore recover much of their read share on PacBio.

### 2.5 Consequences for interpretation

1. Viral recall on both datasets is measured at **~675× (Illumina) and ~65×
   (PacBio) less viral read depth than intended** — 3,336 viral reads per
   Illumina sample and 896–997 per PacBio sample. Every viral figure in this
   analysis is therefore a lower bound on what these pipelines would achieve on
   the community as designed, and the viral results are provisional pending
   re-simulation.
2. **PacBio's higher viral recall is partly a composition effect, not only a
   classifier or read-length-accuracy effect.** The same community delivers
   ~10× more viral read share on PacBio. Platform differences in viral recall
   must not be attributed to long-read classification accuracy alone.
3. The two datasets are approximately **matched in bases, not in reads**:
   16.6 M Illumina read pairs ≈ 5.0 Gb versus 0.45–0.52 M PacBio reads of
   comparable total length. The 33-fold read-count difference is therefore not
   a depth imbalance in the usual sense.
4. This diagnosis is inferred from the outputs, because no CAMISIM
   configuration, seed, genome manifest or log was supplied. Confirming it
   against the actual configuration remains an open item
   (`missing_information_and_questions.md`, A1).

Bacteria and archaea are present at close to the designed 88% of reads (93.8%
Illumina, 98.8% PacBio of total reads), as expected, since the bacterial block
dominates both the copy-number and the length-weighted scales.

---

## 3. File-level anomalies found

### 3.1 PacBio mag-nf read-level reports are the same file three times — EXCLUDED

`Pacbio/mag-nf/sample_{0,1,2}/sample_708_pacbio_{0,1,2}_reads_kraken2_report.txt`
are **byte-identical to each other** (MD5 `70de11486f26…`) and byte-identical to
`Pacbio/TaxProfiler/sample_0/sample_708_0_PacBio_null_Viral_ref.kraken2.kraken2.report.txt`.

These files cannot represent three different PacBio samples. They were excluded
from scoring. PacBio mag-nf was instead scored on its **contig-level** Kraken2
reports, which do differ between samples (98,803 / 100,412 / 101,568 bytes).
This makes the mag-nf Illumina (read-level) versus PacBio (contig-level)
comparison non-equivalent in kind; it is flagged wherever reported.

### 3.2 Illumina TaxTriage and TaxProfiler are the same Kraken2 run

For all three Illumina samples, `Illumina/TaxTriage/sample_N/Illumina_708_N.kraken2.report.txt`
is byte-identical to `Illumina/TaxProfiler/sample_N/sample_708_N_null_Viral_ref.kraken2.kraken2.report.txt`.
Their species-level metrics are therefore identical by construction, not by
independent agreement. The two PacBio runs *are* distinct files and give
slightly different numbers.

### 3.3 Other byte-identical pairs (benign)

- `complete.krakenreport.merged.hierarchy.csv` ≡ `*.kraken2.report.txt.single.csv`
  in all six TaxTriage sample directories (a pipeline alias).
- `Illumina/TaxProfiler/databases.csv` ≡ `Pacbio/TaxProfiler/databases.csv`.
- `Illumina/KrakenMetaReads-nf/metagenome_hybrid_workflow.nf` ≡ the PacBio copy.

### 3.4 Database naming is misleading

`databases.csv` labels the TaxProfiler database `Viral_ref`, but its path is
`.../databases/k2_standard_20240605` — the Kraken2 **Standard** database
(2024-06-05 build), not a viral-only database. This is consistent with
TaxProfiler recovering all 490 bacterial/archaeal members. The same standard
database is used by TaxTriage, mag-nf and KrakenMetaReads-nf. CLARK, GOTTCHA2
and sourmash use different, tool-specific databases (see
`_superseded/tables/Table2_pipelines_versions_databases.csv`).

### 3.5 sourmash uses a different taxonomy namespace

sourmash prokaryote results are GTDB (reps RS226) lineage strings
(`s__Prevotella sp934191715`) with **no NCBI taxids**; viral results come from
an NCBI-viruses 2025.01 sketch database. Its results can therefore only be
matched to the ground truth by normalised name. Where a GTDB species name has no
NCBI counterpart, a non-match may reflect naming systems rather than
non-detection. All sourmash numbers are flagged `match_mode = name` and are not
directly comparable with the taxid-matched pipelines.

### 3.6 No empty, truncated or unreadable files

Every one of the 302 files was opened successfully. The smallest result files
(GOTTCHA2 `.out`, 349–366 bytes) are complete SLURM logs, not truncated data.

---

## 4. Output-format heterogeneity

| Pipeline | Format | Ranks reported | Taxids | Read counts |
|---|---|---|---|---|
| TaxProfiler | Kraken2 report (+ Bracken TSV, Illumina only) | full lineage, U→S2 | yes | yes |
| TaxTriage | Kraken2 report | full lineage | yes | yes |
| mag-nf | Kraken2 report (reads or contigs) | full lineage | yes | yes |
| KrakenMetaReads-nf | Kraken2 report over assembled contigs; Illumina also a merged two-assembler CSV | full lineage | yes | contig-derived counts |
| GOTTCHA2 | own TSV with an explicit `LEVEL` column | superkingdom→species | yes | yes |
| CLARK | abundance CSV with a semicolon lineage string | species only | yes | yes |
| sourmash | `summarized.csv` with GTDB/NCBI lineage strings | superkingdom→species | **no** | weighted fraction only |

Because ranks are mixed, the ground-truth taxid was matched against **any**
species-or-below node (S, S1–S5, species, strain, serotype, no-rank), so that
the 67 phage members (rank S1), the two influenza serotype nodes and the four
no-rank nodes are scorable. This is documented in
`04_name_mapping_and_manual_review.xlsx`, sheet `matching_rules`.

---

## 5. Analysis limits that follow from the audit

1. **No true negatives.** The negative space is each tool's own reference
   database, which differs per tool and is not supplied. Specificity, accuracy,
   balanced accuracy and MCC are therefore **not computed**. Only TP, FP, FN,
   recall, precision, F1, FDR and FNR are reported.
2. **Precision is database-size dependent.** A tool searching the Kraken2
   Standard database can report thousands of taxa absent from a 708-member
   community; a tool with a smaller database cannot. Precision comparisons
   across tools with different databases measure database breadth as much as
   classifier stringency, and are reported alongside a read-threshold sweep
   (≥1, ≥10, ≥100 reads) rather than as a single number.
3. **Illumina and PacBio samples are index-matched but not proven to be paired.**
   `sample_0/1/2` correspond by name across platforms and derive from the same
   708-member design, but no seed or manifest confirms read-level pairing. Paired
   contrasts are reported as sample-wise differences with bootstrap intervals;
   no paired significance test is claimed (n = 3).
4. **Genus- and family-level scoring is not performed**: the ground-truth
   document supplies no genus or family column, and constructing one from tool
   output would make the ground truth tool-dependent.
5. **Pipeline versions are read from run scripts and logs.** Where a version is
   not recorded in any supplied file, it is marked as such rather than inferred.

---

## 6. Files excluded from scoring, by reason (summary)

| Reason | n files |
|---|---|
| Derived summary; primary table scored instead | 21 |
| Visualisation input (Krona), redundant with the Kraken2 report | 18 |
| Multi-sample combined report; single-sample report scored instead | 15 |
| No NCBI taxids; `summarized.csv` scored instead (sourmash `.kreport.txt`) | 12 |
| Redundant with `summarized.csv` (sourmash `.human.txt`) | 12 |
| Truncated top-N view of the abundance table | 6 |
| GOTTCHA2 `.full.tsv` superseded by the filtered `.tsv` | 6 |
| GOTTCHA2 `.lineage.tsv` redundant with the scored `.tsv` | 6 |
| No read counts; Kraken2 report scored instead (TaxTriage `.single.csv`) | 6 |
| Contig-level, carries no taxon assignment | 6 |
| Byte-identical to the sibling `.single.csv` | 6 |
| Contig-level; the metaFlye Kraken2 report scored instead | 3 |
| 3-line text summary; `merged_report.csv` scored instead | 3 |
| Bracken stream taken from `.bracken.tsv` instead | 3 |
| Read-level assignment file; aggregated form scored instead | 3 |
| Logs, scripts, configs, sample sheets, task brief | 95 |
| Source documents used as ground truth / background | 3 |

Per-file reasons are in `01_file_inventory.xlsx`, column `reason_if_excluded`.
