#!/bin/bash
#SBATCH --job-name=k2_pacbio_sam1
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=k2_pacbio_sam1_%j.out
#SBATCH --error=k2_pacbio_sam1_%j.err

# ==============================================================================
# Taxonomic classification for the PacBio HiFi mag run.
#
# WHY: nf-core/mag (as configured here) only ASSEMBLES (Flye) — it has no
# Kraken2 DB, so it produces NO species classification. This script adds the
# classification step, WITHOUT re-assembling:
#   1) Kraken2 on the Flye contigs that mag already produced   (contig-level)
#   2) Kraken2 on the raw PacBio HiFi reads                    (read-level, more
#                                                               sensitive to rare taxa)
# Both use the same Standard Kraken2 DB as the rest of the benchmark.
#
# Submit:  sbatch run_kraken2_classify_pacbio.sh
# ==============================================================================

set -Eeuo pipefail
timestamp() { date +"%Y-%m-%d %H:%M:%S"; }
log()  { echo "[$(timestamp)] $*"; }
die()  { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
# Guard against a stale/inherited SLURM_SUBMIT_DIR (e.g. an OOD dashboard
# session env var pointing at /var/www/ood/apps/sys/dashboard) that does not
# exist on this filesystem. Fall back to $PWD in that case instead of failing
# the cd outright.
[[ -d "$SUBMIT_DIR" ]] || { echo "WARNING: SLURM_SUBMIT_DIR=$SUBMIT_DIR does not exist, falling back to \$PWD" >&2; SUBMIT_DIR="$PWD"; }
cd "$SUBMIT_DIR"

# ------------------------------------------------------------------
# Config (override at submission with --export=ALL,VAR=...) Edited:41,42,44
# ------------------------------------------------------------------
SAMPLE_NAME="${SAMPLE_NAME:-sample_708_pacbio_2}"
MAG_OUTDIR="${MAG_OUTDIR:-results_pacbio_hifi_sample_2}"
# Raw PacBio HiFi reads (read-level classification).
RAW_READS="${RAW_READS:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/pacbio-708-genomes/pacbio_hifi_reads_sample_2.fastq.gz}"
KRAKEN2_DB="${KRAKEN2_DB:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/k2_standard_20240605}"
KRAKEN2_CONTAINER="${KRAKEN2_CONTAINER:-docker://quay.io/biocontainers/kraken2:2.1.2--pl5262h7d875b9_0}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/tmp}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

OUTDIR="${SUBMIT_DIR}/kraken2_classification_${SAMPLE_NAME}"
mkdir -p "$OUTDIR"

# ------------------------------------------------------------------
# Locate the Flye contigs mag already produced (published under Assembly/).
# nf-core/mag writes e.g. Assembly/FLYE/FLYE-<sample>.assembly.fasta.gz
# ------------------------------------------------------------------
command -v apptainer >/dev/null 2>&1 || die "apptainer not found"
[[ -d "$KRAKEN2_DB" ]]  || die "Kraken2 DB not found: $KRAKEN2_DB"
[[ -f "$RAW_READS" ]]   || die "Raw reads not found: $RAW_READS"
[[ -d "$MAG_OUTDIR" ]]  || die "mag outdir not found: $MAG_OUTDIR (set MAG_OUTDIR=...)"

CONTIGS="$(find "$MAG_OUTDIR" -type f \( -iname '*FLYE*assembly.fasta.gz' -o -iname '*assembly.fasta.gz' \) 2>/dev/null | head -1)"
[[ -n "$CONTIGS" && -s "$CONTIGS" ]] || die "Could not find Flye contigs under $MAG_OUTDIR (looked for *assembly.fasta.gz). Set CONTIGS=... explicitly."
log "Using mag Flye contigs: $CONTIGS"

CONTIG_REPORT="${OUTDIR}/${SAMPLE_NAME}_contigs_kraken2_report.txt"
CONTIG_OUT="${OUTDIR}/${SAMPLE_NAME}_contigs_classification.txt"
READ_REPORT="${OUTDIR}/${SAMPLE_NAME}_reads_kraken2_report.txt"
READ_OUT="${OUTDIR}/${SAMPLE_NAME}_reads_classification.txt"

log "=========================================="
log "Sample     : $SAMPLE_NAME"
log "Contigs    : $CONTIGS"
log "Raw reads  : $RAW_READS"
log "Kraken2 DB : $KRAKEN2_DB"
log "Output     : $OUTDIR"
log "=========================================="

# ------------------------------------------------------------------
# 1) Kraken2 on the mag Flye contigs (gzipped fasta -> --gzip-compressed)
# ------------------------------------------------------------------
log "Step 1 - Kraken2 on Flye contigs"
apptainer exec --no-mount /lscratch --bind /scratch:/scratch "$KRAKEN2_CONTAINER" \
    kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" --use-names --gzip-compressed \
            --output "$CONTIG_OUT" --report "$CONTIG_REPORT" "$CONTIGS" \
    || die "Kraken2 on contigs failed"
[[ -s "$CONTIG_REPORT" ]] || die "empty contig Kraken2 report"
log "  -> $CONTIG_REPORT"

# ------------------------------------------------------------------
# 2) Kraken2 on the raw PacBio HiFi reads (read-level, most sensitive)
# ------------------------------------------------------------------
log "Step 2 - Kraken2 on raw reads"
apptainer exec --no-mount /lscratch --bind /scratch:/scratch "$KRAKEN2_CONTAINER" \
    kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" --use-names --gzip-compressed \
            --output "$READ_OUT" --report "$READ_REPORT" "$RAW_READS" \
    || die "Kraken2 on reads failed"
[[ -s "$READ_REPORT" ]] || die "empty reads Kraken2 report"
log "  -> $READ_REPORT"

# ------------------------------------------------------------------
# Summary: species-level lines (rank code 'S') from each report
# ------------------------------------------------------------------
n_contig_sp=$(awk -F'\t' '$4=="S"' "$CONTIG_REPORT" | wc -l)
n_read_sp=$(awk -F'\t' '$4=="S"' "$READ_REPORT" | wc -l)
log "=========================================="
log "Done - $SAMPLE_NAME"
log "  contig-level species (rank S): $n_contig_sp   -> $CONTIG_REPORT"
log "  read-level   species (rank S): $n_read_sp   -> $READ_REPORT"
log "Top read-level species:"
awk -F'\t' '$4=="S"{print $2"\t"$6}' "$READ_REPORT" | sort -t$'\t' -k1 -nr | head -10 | sed 's/^/    /'
log "End: $(timestamp)"
