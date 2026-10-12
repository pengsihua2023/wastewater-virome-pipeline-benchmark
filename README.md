# Seven-pipeline viral-detection benchmark on simulated wastewater (Illumina vs PacBio HiFi)

Code and data supporting: *[manuscript title]*, submitted to Water Research.
Archive DOI: **10.5281/zenodo.XXXXXXX** (release `v1.0-submission`). GitHub: pengsihua2023/wastewater-virome-pipeline-benchmark

## What this is
A benchmark of seven metagenomic pipelines (nf-core/taxprofiler, nf-core/taxtriage, nf-core/mag, KrakenMetaReads-nf, GOTTCHA2, CLARK, sourmash) on simulated Illumina (CAMISIM) and PacBio HiFi (PBSIM3 + ccs) wastewater data. The 708-member community (151 non-phage viruses = ground truth, 67 phages, 490 bacteria/archaea) is derived from the public wastewater metatranscriptome HTP-2024-01-15 (Los Angeles, Grimm et al.). No human-subject data were generated or analysed.

## Layout
| Path | Content |
|---|---|
| `code/analysis/` | numbered analysis scripts (run in numeric order; 31–33 produce the threshold, community and revision tables) |
| `code/simulation/` | read-simulation code: `simulate_pacbio.py`, `run_pacbio_sample.sh`, `run_art*.sh`, `merge_and_verify.sh`; `versions/` holds all historical versions of the PacBio script (see its README for the version that produced the data) |
| `data/regen/` | simulation design: genome lengths, per-sample read plans, ART/PBSIM task tables, per-member read counts, md5 manifests |
| `data/genome_locations.tsv`, `data/distribution_{0,1,2}.txt` | genome list and per-sample abundance files used for the simulation |
| `data/analysis_708/interim/` | parsed detections and metrics; `interim/coverage/` holds the sources of Supplementary Tables S1–S14 |
| `data/analysis_708/figures_*` | manuscript figures |
| `data/Supplementary_Tables_coverage.xlsx` | Supplementary Tables S1–S14 (README sheet describes each) |
| `data/Illumina/`, `data/Pacbio/` | raw outputs of the 42 pipeline runs (**not in this repository; download from the Zenodo archive and unpack into `data/`**) |
| `docs/` | regeneration guide and the verification report of the simulated data |

All scripts take `data/` as their root; run them from anywhere, e.g. `python code/analysis/31_per_virus_threshold_detection.py`.

## Reproducing the results
1. `python code/analysis/02_parse_pipeline_outputs.py` → parse the 42 runs (needs `data/Illumina/`, `data/Pacbio/` from Zenodo)
2. `code/analysis/14_viral_only_scoring.py`, `15_genus_level_sensitivity.py` → viral-only scoring (taxid match at species-or-below; thresholds ≥1/≥10/≥100 clade reads)
3. `code/analysis/31`, `32`, `33`, then `28_floor_recall_and_supp_tables.py` → per-virus detection, full community table, bootstrap, supplementary workbook (the large parsed table `data/analysis_708/interim/viral_detections_long.tsv` is provided on Zenodo, not in this repository)
4. `code/analysis/17_platform_figures.py`, `18_combined_figures.py` → figures

Requirements: Python ≥ 3.10, pandas, numpy, openpyxl, matplotlib, python-docx (see `requirements.txt`). Random seed for bootstrap: 20261012.

## Reference genomes (`genomes/`, not archived, 2.3 GB)
The `path` column of `data/regen/*.tsv` and `data/genome_locations.tsv` is relative to a `genomes/` directory that holds the 708 community genomes (`genomes/Viruses/`, `genomes/Bacteriophage/`, `genomes/Bacteria_Archaea/`, 151 + 67 + 490 FASTA files, NCBI RefSeq/GenBank assemblies named `<taxid>_<accession>_genomic.fna`). Download them from NCBI by accession and place them there; `data/regen/genomes_md5.txt` lists the md5 of every file used so the set can be verified (`cd <parent of genomes/> && md5sum -c <repo>/data/regen/genomes_md5.txt`). Scripts that read genome files (e.g. `32_full_community_table.py`) need this directory; scripts that only read the TSV tables do not.

## Regenerating the simulated reads (not archived, ~29 GB)
Use `data/regen/plan_*.tsv` with `code/simulation/run_art.sh` (Illumina) and `code/simulation/simulate_pacbio.py` (PacBio); the shell scripts use paths relative to `data/`, so run them from inside `data/`. Expected per-member read counts are in `data/regen/per_member_s*.tsv`. md5 checksums of the original FASTQ files (6 Illumina, 3 PacBio HiFi) are in `data/regen/fastq_md5.txt`; regenerated reads are expected to match in read counts, not necessarily byte for byte. See `docs/guide-line-regenerat-data.md`.

## Software and databases
Versions of all pipelines, Kraken2/Bracken, CLARK, GOTTCHA2 (locally patched for PacBio), sourmash and databases: [to add — see Methods §2.2]; ccs version: [to add].

## Notes
- Strong simplifications apply: three simulated replicates share the same genomes and abundance design; see the manuscript limitations.
- Some very large files (CLARK per-read assignments, long-format detection tables) are provided as separate Zenodo files [to confirm].

## License and citation
Code: MIT. Data and tables: CC BY 4.0. Cite the manuscript and the Zenodo DOI (`CITATION.cff`).

## Contact
Sihua Peng, College of Public Health, University of Georgia.
