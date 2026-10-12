# Pipeline run scripts and configurations

Submission scripts (SLURM), Nextflow configurations and sample sheets used to run the seven pipelines on the six simulated datasets
(three Illumina, three PacBio HiFi). Paths were generalised for publication: the HPC user name is replaced by `$USER`, and
site-specific directories (`/scratch/...`, database locations) must be adapted to your system. The scripts are provided as a record of the
exact parameters used; they are not packaged as a portable workflow. Comments of the form "Edited: ..." are the author's notes on which lines were changed between runs.

| Directory | Pipeline | Files |
|---|---|---|
| `Illumina/TaxProfiler/`, `PacBio/TaxProfiler/` | nf-core/taxprofiler v1.2.0 (Kraken2) | run/submit scripts, `nextflow_*.config`, sample sheets, abundance calculation scripts (`calculate_abundance*_en.py`, `batch_calculate_abundance*_en.sh`) |
| `Illumina/TaxTriage/`, `PacBio/TaxTriage/` | nf-core/taxtriage v1.0 | submit scripts, sample sheets |
| `Illumina/mag-nf/`, `PacBio/mag-nf/` | nf-core/mag v3.1.0 | run scripts, configs, sample sheets |
| `Illumina/KrakenMetaReads/`, `PacBio/KrakenMetaReads/` | KrakenMetaReads-nf (in-house Nextflow DSL2) | `metagenome_hybrid_workflow.nf` (identical in both directories), configs, run scripts, sample sheets |
| `Illumina/GOTTCHA2/`, `PacBio/GOTTCHA2/` | GOTTCHA2 | run scripts (PacBio: locally patched binary, see the paper) |
| `Illumina/CLARK/`, `PacBio/CLARK/` | CLARK | SLURM classification scripts |
| `Illumina/sourmash/`, `PacBio/sourmash/` | sourmash 4.9.4 | run scripts |

## KrakenMetaReads-nf: version used here versus the public repository

The in-house workflow is also published at https://github.com/pengsihua2023/KrakenMetaReads-nf (workflow version 3.0.0).
The `metagenome_hybrid_workflow.nf` in this directory is the version that was **actually used for the benchmark** and differs from the
public repository version: it renames reads to unique identifiers before Flye assembly, calls metaSPAdes with explicit k-mer sizes
(`-k 21,33,55,77,99,127`), and adds `--genome-size` and `--min-overlap` options and a modified output (publishDir) definition for the Flye step.
For the benchmark results, use the file archived here.
