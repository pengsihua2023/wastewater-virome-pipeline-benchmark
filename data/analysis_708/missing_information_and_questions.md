# Information required before a submission-ready version

Ordered by how much each item affects the conclusions.

---

## A. Blocking for the viral conclusions

### A1. RESOLVED as to intent — the datasets need re-simulating

**Resolved:** the 88 / 8 / 4 split is intended as a share of **reads**. This is
confirmed by the study authors and is consistent with the design document's
source column, which holds Kraken2 read counts. Among the 708 selected taxa in
sample HTP-2024-01-15 the real read split was 5.6% bacteria/archaea, 1.0% phage
and 93.4% virus, so a 4% viral read target is a deliberate and conservative
downscaling of the source metatranscriptome, not an inflation of it.

**Note also** that the document specifies the viral target in reads twice, and
inconsistently: the `Simulation design reads` column implies a 5.6 / 1.0 / 93.4%
read split (it is the real sample's clade reads), while `Relative abundance in
simulation (%)` renormalises within each block to 88 / 8 / 4. Please confirm
which is authoritative and remove the other, or label them clearly.

**Outstanding problem:** CAMISIM consumes its abundance file as genome copies,
so the 88 / 8 / 4 input produced a viral **read** share of 0.020% (Illumina) and
0.18–0.22% (PacBio) — missing the 12% target by ~675× and ~65× respectively, and
missing the per-member design read counts by 10^4 or more (*Mamastrovirus 1*
559,347 specified vs 606 delivered). Two independent tests establish this:
per-domain read yield per unit abundance matches the genome-length prediction to
within ~20% on Illumina (`scripts/06_abundance_scale_diagnosis.py`), and
length-weighting the two candidate input columns predicts 2.84% versus 0.017%
viral read share against 0.017% measured, identifying the 88 / 8 / 4 column as
the input (`scripts/07_designed_vs_observed_reads.py`).

**Action required — re-simulation.** The copy abundances that would deliver the
intended read shares are in `interim/abundance_correction_required.tsv`:

| Platform | Bacteria/archaea | Bacteriophage | Virus |
|---|---|---|---|
| Current design | 88.0% | 8.0% | 4.0% |
| Required (Illumina) | 0.93% | 71.7% | 27.3% |
| Required (PacBio) | 7.7% | 77.0% | 15.3% |

Two design decisions follow:

1. The required abundance file **differs between platforms** (read-length capping
   changes the yield per copy for small genomes), so one file cannot give both
   platforms a 12% viral read share. Either accept different viral read shares
   per platform, or generate a platform-specific abundance file and state so.
2. The general fix is to set each member's copy abundance to
   (target read share) ÷ (genome length). This needs a genome-length column,
   which the design document does not currently carry — please add it.

Still needed to confirm the diagnosis directly rather than by inference from
pipeline outputs: the CAMISIM configuration file(s), the abundance/distribution
file actually passed to CAMISIM, the genome manifest (accession per member, and
whether any download failed), the random seed(s) and the simulation logs.

**Until re-simulation, every viral figure in this analysis is a lower bound**,
measured at roughly two orders of magnitude less viral read depth than intended.

### A2. Were the Illumina and PacBio `sample_N` datasets generated from the same draws?

If they were generated from the same seed and the same per-genome abundances,
platform contrasts can be treated as paired and analysed accordingly. Nothing in
the supplied files establishes this, so the analysis reports unpaired
sample-wise differences with bootstrap intervals and claims no paired test.

### A3. Library sizes are base-matched — confirm this was intentional

Illumina samples contain ~16.6 M read pairs (≈5.0 Gb) and PacBio samples
0.45–0.52 M long reads of comparable total length, so the datasets appear to be
matched in **bases** rather than reads. Please confirm. If so, the 33-fold
read-count difference is not a depth imbalance and should be described as
base-matched in the Methods.

Separately: for the same base-matched community, PacBio carries ~10× the viral
read share of Illumina (Section A1). Whether the platform comparison should be
reported at equal viral read depth rather than at equal total depth is a design
decision for the authors.

---

## B. Blocking for the Methods section

### B1. PacBio nf-core/mag read-level reports are the same file three times

`Pacbio/mag-nf/sample_{0,1,2}/sample_708_pacbio_{0,1,2}_reads_kraken2_report.txt`
are byte-identical to each other and to the PacBio taxprofiler report for
`sample_0` (MD5 `70de11486f26…`). Was the read-level Kraken2 step re-run per
sample, or were these files copied? If per-sample reports exist, they should
replace the contig-level evidence currently used for PacBio mag-nf.

### B2. Illumina TaxTriage and TaxProfiler share one Kraken2 report

For all three Illumina samples the two pipelines' Kraken2 reports are
byte-identical. Were these genuinely separate runs whose outputs coincide
byte-for-byte (unlikely), or was one report copied into the other pipeline's
output directory? As it stands, Illumina TaxTriage is not an independent
evaluation.

### B3. Missing software versions

Not recorded in any supplied file, currently unstated in the manuscript:
- CAMISIM version and the read simulators it invoked (ART? PBSIM? which model?)
- GOTTCHA2 version and signature-database build date
- CLARK version and CLARK_DB build date
- exact nf-core/taxtriage revision (the run uses a local checkout, not `-r`)
- Nextflow, Java and container-runtime versions actually used
- Kraken2 Standard database: confirmed as `k2_standard_20240605`; the CLARK,
  GOTTCHA2 and sourmash database build dates are not recorded

### B4. Database labelling

`databases.csv` names the taxprofiler database `Viral_ref`, but its path is the
Kraken2 **Standard** database. Confirm that this is intended, since the label
would otherwise mislead readers about the pipeline's scope.

### B5. GOTTCHA2 PacBio patch

The PacBio run script applies a local patch to `gottcha2.py` to add `map-hifi`
to the `-xm` choices. Confirm the patch content and whether upstream GOTTCHA2
supports HiFi input, since running an unreleased modification affects
reproducibility and the interpretation of GOTTCHA2's PacBio results.

### B6. sourmash Illumina input — checked, no action needed

The `query_filename` field of the sourmash `summarized.csv` files lists only
`sample_708_N_02.fq.gz`, which initially suggested that only R2 had been
sketched. Inspection of `Illumina/sourmash/run_sourmash_Illumina.sh` (lines
88–92) confirms that both R1 and R2 are passed to `sourmash sketch dna`; the
field simply records the last input. No correction is required.

---

## C. Needed for the manuscript front and back matter

- Author list, affiliations and CRediT contribution statement.
- Funding sources and grant numbers.
- Conflict-of-interest declarations.
- Whether the simulated datasets will be deposited (SRA/Zenodo) and, if so,
  their accessions.
- A public repository URL and archived DOI for the analysis code and for the
  in-house KrakenMetaReads-nf workflow.
- Target journal, so reference style and figure specifications can be matched.

---

## D. References requiring verification

Cited by name and version in the Methods but not verifiable from the project
files; each must be checked against the publication record before submission:

- nf-core/taxprofiler, nf-core/taxtriage and nf-core/mag pipeline papers
- metaSPAdes, MEGAHIT, fastp, minimap2, samtools primary references
- GTDB release RS226 reference
- NCBI viruses 2025.01 sourmash database reference

No reference has been invented; anything not confirmed from the supplied
documents is marked as pending in the manuscript reference list.

---

## E. Scientific questions the authors may wish to answer

1. Should the 13 manually added sentinel viruses be scored separately from the
   138 automatically selected ones? They were assigned abundances by expert
   judgement rather than measured from the sample, so they test a different
   thing.
2. Is a genus- or family-level evaluation wanted? It was not performed because
   the ground-truth document carries no genus or family column; supplying one
   would enable it.
3. Should CLARK's γ ≥ 0.2 filter be the primary CLARK result rather than a
   sensitivity analysis? The current choice (unfiltered) matches what the other
   pipelines report but is not what a CLARK user would deploy.
4. Should a minimum-read threshold be adopted as the primary operating point?
   The ≥ 10 and ≥ 100 read results are already computed and change the precision
   ranking substantially.
