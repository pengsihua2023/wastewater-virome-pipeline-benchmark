# Seven-pipeline viral-detection benchmark on simulated wastewater (Illumina vs PacBio HiFi)

Code and data supporting: *[manuscript title]*, submitted to Water Research.
Archive DOI: **10.5281/zenodo.XXXXXXX** (release `v1.0-submission`). GitHub: pengsihua2023/wastewater-virome-pipeline-benchmark

## What this is
A benchmark of seven metagenomic pipelines (nf-core/taxprofiler, nf-core/taxtriage, nf-core/mag, KrakenMetaReads-nf, GOTTCHA2, CLARK, sourmash) on simulated Illumina (CAMISIM) and PacBio HiFi (PBSIM3 + ccs) wastewater data. The 708-member community (151 non-phage viruses = ground truth, 67 phages, 490 bacteria/archaea) is derived from the public wastewater metatranscriptome HTP-2024-01-15 (Los Angeles, Grimm et al.). No human-subject data were generated or analysed.

## Layout
| Path | Content |
|---|---|
| `regen/` | simulation design: genome lengths, per-sample read plans, ART/PBSIM task tables, per-member read counts, helper scripts |
| `simulate_pacbio.py`, `run_pacbio_sample.sh`, `simulation-code-versions/` | PacBio simulation code (see its README for the version used to produce the data) |
| `Illumina/<pipeline>/sample_{0,1,2}/`, `Pacbio/<pipeline>/sample_{0,1,2}/` | raw outputs of the 42 pipeline runs |
| `analysis_708/scripts/` | numbered analysis scripts (run in numeric order; 31–33 produce the threshold, community and revision tables) |
| `analysis_708/interim/` | parsed detections and metrics; `interim/coverage/` holds the sources of Supplementary Tables S1–S14 |
| `analysis_708/figures_*` | manuscript figures |
| `Supplementary_Tables_coverage.xlsx` | Supplementary Tables S1–S14 (README sheet describes each) |

## Reproducing the results
1. `python analysis_708/scripts/02_parse_pipeline_outputs.py` → parse the 42 runs
2. `…/14_viral_only_scoring.py`, `15_genus_level_sensitivity.py` → viral-only scoring (taxid match at species-or-below; thresholds ≥1/≥10/≥100 clade reads)
3. `…/31`, `32`, `33`, then `28_floor_recall_and_supp_tables.py` → per-virus detection, full community table, bootstrap, supplementary workbook
4. `…/17_platform_figures.py`, `18_combined_figures.py` → figures

Requirements: Python ≥ 3.10, pandas, numpy, openpyxl, matplotlib, python-docx (see `requirements.txt`). Random seed for bootstrap: 20261012.

## Reference genomes (`genomes/`, not archived, 2.3 GB)
The `path` column of `regen/*.tsv` and `genome_locations.tsv` is relative to a `genomes/` directory that holds the 708 community genomes (`genomes/Viruses/`, `genomes/Bacteriophage/`, `genomes/Bacteria_Archaea/`, 151 + 67 + 490 FASTA files, NCBI RefSeq/GenBank assemblies named `<taxid>_<accession>_genomic.fna`). Download them from NCBI by accession and place them there; `regen/genomes_md5.txt` lists the md5 of every file used so the set can be verified (`cd <parent of genomes/> && md5sum -c regen/genomes_md5.txt`). Scripts that read genome files (e.g. `32_full_community_table.py`) need this directory; scripts that only read the TSV tables do not.

## Regenerating the simulated reads (not archived, ~29 GB)
Use `regen/plan_*.tsv` with `run_art.sh` (Illumina) and `simulate_pacbio.py` (PacBio); expected per-member read counts are in `regen/per_member_s*.tsv`. md5 checksums of the original FASTQ files (6 Illumina, 3 PacBio HiFi) are in `regen/fastq_md5.txt`; regenerated reads are expected to match in read counts (`regen/per_member_s*.tsv`), not necessarily byte for byte.

## Software and databases
Versions of all pipelines, Kraken2/Bracken, CLARK, GOTTCHA2 (locally patched for PacBio), sourmash and databases: [to add — see Methods §2.2]; ccs version: [to add].

## Notes
- Strong simplifications apply: three simulated replicates share the same genomes and abundance design; see the manuscript limitations.
- Some very large files (CLARK per-read assignments, long-format detection tables) are provided as separate Zenodo files [to confirm].

## License and citation
Code: MIT. Data and tables: CC BY 4.0. Cite the manuscript and the Zenodo DOI (`CITATION.cff`).

## Contact
Sihua Peng, College of Public Health, University of Georgia.
