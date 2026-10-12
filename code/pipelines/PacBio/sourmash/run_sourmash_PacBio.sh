#!/bin/bash

#SBATCH --job-name=Sourmash_sam1-PacBio
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=5-04:00:00
#SBATCH --output=Sourmash_sam1-PacBio_%j.out
#SBATCH --error=Sourmash_sam1-PacBio_%j.err

# ============================================================================
# PacBio HiFi long reads -> metaFlye assembly -> sourmash qualitative virome classification
#   Workflow: raw HiFi reads --(metaFlye --pacbio-hifi --meta)--> contig
#          --(sourmash sketch/gather/tax)--> which viruses/phages are detected
#   Note: long-read assembly follows "assemble first, then classify", avoiding the high error rate of raw reads
#         The result is [qualitative] (which viruses were detected), not quantitative abundance (contigs carry no abundance information)
# ============================================================================

set -euo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1

# ----------------------------------------------------------------------------
# 0. parameters -- change to your actual paths/samples
## Edited: 27, 29, 45
# ----------------------------------------------------------------------------
SAMPLE="pacbio_hifi_reads_sample_1"                    # sample name (output prefix)
# Input: raw PacBio HiFi reads (single file, long reads have no R1/R2)
READS_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/pacbio-708-genomes"
READS="${READS_DIR}/${SAMPLE}.fastq.gz"   # <-- change to your actual file name

# Database
DB_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/sourmash_db"
VIRAL_DB="${DB_DIR}/ncbi-viruses-2025.01.k31.zip"
VIRAL_LIN="${DB_DIR}/ncbi-viruses.2025.01.lineages.csv"
TAX_SQLDB="${DB_DIR}/ncbi-viruses.2025.01.taxonomy.sqldb"

# sourmash parameters (must match the database: k=31, scaled=50)
# HiFi error rate <1%, close to Illumina, k=31 is suitable
KSIZE=31
SCALED=50
THRESHOLD_BP=100

# Output directory (placed under scratch, writable by the container)
OUTBASE="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/sourmash/Pacbio-708/results_pacbio_708_1"
OUTDIR="${OUTBASE}/${SAMPLE}"
ASMDIR="${OUTDIR}/metaflye"               # metaFlye Output directory
mkdir -p "${OUTDIR}"

THREADS=${SLURM_CPUS_PER_TASK:-16}

# ----------------------------------------------------------------------------
# 1. Load the software (Apptainer container)
# ----------------------------------------------------------------------------
SIF_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/containers"
SOURMASH_SIF="${SIF_DIR}/sourmash_4.9.4.sif"
FLYE_SIF="${SIF_DIR}/flye_2.9.5.sif"      # <-- must be pulled first, see the note below

# Wrap as a function so that the script can call sourmash / flye
sourmash() { apptainer exec --bind /scratch "${SOURMASH_SIF}" sourmash "$@"; }
flye()     { apptainer exec --bind /scratch "${FLYE_SIF}"     flye     "$@"; }

# if apptainer on this node needs module load first, uncomment::
# module load apptainer

# ---- How to pull the flye container (one time, run on a node where apptainer can be used) ----
#   cd ${SIF_DIR}
#   apptainer pull flye_2.9.5.sif docker://quay.io/biocontainers/flye:2.9.5--py310h275bdba_3
#   apptainer exec flye_2.9.5.sif flye --version   # Verify

echo "==== $(date) | sample: ${SAMPLE} (PacBio HiFi) ===="

# Check that the input file exists (catches path errors early)
if [[ ! -s "${READS}" ]]; then
    echo "ERROR: Input reads do not exist or are empty: ${READS}" >&2
    echo "       Please check the READS_DIR and SAMPLE naming" >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# 1.5. Preprocessing: remove duplicated read ID
#      CAMISIM simulated data may contain duplicate read IDs (e.g., S1/1/ccs appears repeatedly)
#      Flye will refuse to process input with duplicated IDs, so the IDs must be made unique first
# ----------------------------------------------------------------------------
READS_DEDUP="${ASMDIR}/reads_dedup.fastq.gz"
if [[ ! -s "${READS_DEDUP}" ]]; then
    echo "---- preprocessing: make read IDs unique ----"
    mkdir -p "${ASMDIR}"
    zcat "${READS}" | \
        awk 'BEGIN {count=0} \
             NR % 4 == 1 {count++; $0 = "@read_" count; print; next} \
             {print}' | \
        gzip > "${READS_DEDUP}"
    echo "Deduplicated reads: ${READS_DEDUP}"
else
    echo "---- deduplicated reads already exist, skipping ----"
fi

# Verify the deduplication result
DEDUP_COUNT=$(zcat "${READS_DEDUP}" | grep -c "^@read_")
echo "Reads after processing: ${DEDUP_COUNT}"

# ----------------------------------------------------------------------------
# 2. metaFlye assembly (PacBio HiFi mode)
#    --pacbio-hifi : HiFi specific parameter (error rate <1%)
#    --meta        : metagenome / uneven coverage mode
#    --plasmids option has been deprecated and merged; newer versions automatically try to capture small sequences (viruses/plasmids)
# ----------------------------------------------------------------------------
CONTIGS="${ASMDIR}/assembly.fasta"
if [[ ! -s "${CONTIGS}" ]]; then
    echo "---- metaFlye assembly (--pacbio-hifi --meta) ----"
    flye --version
    flye \
        --pacbio-hifi "${READS_DEDUP}" \
        --meta \
        --threads ${THREADS} \
        --out-dir "${ASMDIR}"
else
    echo "---- assembly already exists, skipping metaFlye ----"
fi

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: metaFlye produced no contigs (${CONTIGS}); reads may be too few or of too poor quality" >&2
    exit 1
fi
echo "Assembly finished: $(grep -c '^>' ${CONTIGS}) sequences: contig"

# ----------------------------------------------------------------------------
# 3. Prepare the taxonomy sqldb (needed only once)
# ----------------------------------------------------------------------------
if [[ ! -s "${TAX_SQLDB}" ]]; then
    echo "---- Building taxonomy sqldb ----"
    sourmash tax prepare --taxonomy-csv "${VIRAL_LIN}" -o "${TAX_SQLDB}" -F sql
else
    echo "---- taxonomy sqldb already exists, skipping ----"
fi

# ----------------------------------------------------------------------------
# 4. SKETCH -- on the contigs (note: noabund, abundance is meaningless for contigs)
# ----------------------------------------------------------------------------
SIG="${OUTDIR}/${SAMPLE}.contigs.sig.zip"
echo "---- sketch (contig) ----"
sourmash sketch dna \
    -p k=${KSIZE},scaled=${SCALED},noabund \
    --name "${SAMPLE}" \
    -o "${SIG}" \
    "${CONTIGS}"
sourmash sig fileinfo "${SIG}"

# ----------------------------------------------------------------------------
# 5. GATHER -- against the viral database
# ----------------------------------------------------------------------------
GATHER_CSV="${OUTDIR}/${SAMPLE}.viral.gather.csv"
echo "---- gather ----"
set +e
sourmash gather "${SIG}" "${VIRAL_DB}" \
    -k ${KSIZE} \
    --scaled ${SCALED} \
    --threshold-bp ${THRESHOLD_BP} \
    --estimate-ani-ci \
    -o "${GATHER_CSV}" \
    --output-unassigned "${OUTDIR}/${SAMPLE}.viral.unassigned.sig.zip"
GATHER_RC=$?
set -e

if [[ ${GATHER_RC} -ne 0 && ! -s "${GATHER_CSV}" ]]; then
    echo "WARNING: sample ${SAMPLE}: no virus detected (gather found no match)"
    echo "==== $(date) | done (no viral hits) ===="
    exit 0
fi

# ----------------------------------------------------------------------------
# 6. TAX METAGENOME -- qualitative virus/phage classification (with full lineage)
#    contig has no abundance, use --no-abundances; the result is "which viruses were detected", not quantitative
# ----------------------------------------------------------------------------
echo "---- tax metagenome (qualitative) ----"
sourmash tax metagenome \
    --gather-csv "${GATHER_CSV}" \
    --taxonomy "${TAX_SQLDB}" \
    --output-base "${OUTDIR}/${SAMPLE}.viral" \
    --output-format csv_summary krona kreport human \
    --no-abundances \
    --rank species

# ----------------------------------------------------------------------------
# 7. Bacteria + archaea (GTDB rs226) -- independent gather
#    The viral database is scaled=50 and the GTDB k31 library is usually scaled=1000; mixing them in one gather would force both to the same value,
#    1000 with a loss of viral sensitivity, so GTDB is gathered separately for the same SIG.
#    scaled=50 sketch can be downsampled to 1000 (FracMinHash nesting), so it is compatible.
# ----------------------------------------------------------------------------
GTDB_DB="${DB_DIR}/gtdb-reps-rs226-k31.dna.zip"                 # <-- confirm the actual zip name/path
GTDB_LIN="${DB_DIR}/gtdb-rs226.lineages.csv"
GTDB_TAX_SQLDB="${DB_DIR}/gtdb-rs226.taxonomy.sqldb"
GTDB_SCALED=1000                                        # GTDB k31 libraries are usually scaled=1000

if [[ ! -s "${GTDB_DB}" ]]; then
    echo "WARNING: GTDB database ${GTDB_DB} not found, skipping bacteria/archaea detection" >&2
else
    if [[ ! -s "${GTDB_TAX_SQLDB}" ]]; then
        echo "---- Building GTDB taxonomy sqldb ----"
        sourmash tax prepare --taxonomy-csv "${GTDB_LIN}" -o "${GTDB_TAX_SQLDB}" -F sql
    fi

    GTDB_GATHER_CSV="${OUTDIR}/${SAMPLE}.prok.gather.csv"
    echo "---- gather (GTDB Bacteria/archaea, scaled=${GTDB_SCALED}) ----"
    set +e
    sourmash gather "${SIG}" "${GTDB_DB}" \
        -k ${KSIZE} \
        --scaled ${GTDB_SCALED} \
        --threshold-bp ${THRESHOLD_BP} \
        --estimate-ani-ci \
        -o "${GTDB_GATHER_CSV}"
    set -e

    if [[ -s "${GTDB_GATHER_CSV}" ]]; then
        echo "---- tax metagenome (GTDB, qualitative) ----"
        sourmash tax metagenome \
            --gather-csv "${GTDB_GATHER_CSV}" \
            --taxonomy "${GTDB_TAX_SQLDB}" \
            --output-base "${OUTDIR}/${SAMPLE}.prok" \
            --output-format csv_summary krona kreport human \
            --no-abundances \
            --rank species
        echo "Bacteria/archaea results: ${OUTDIR}/${SAMPLE}.prok.summarized.csv"
    else
        echo "WARNING: GTDB gather no match (${SAMPLE})"
    fi
fi

echo "==== $(date) | Done ===="
echo "Detected viruses list: ${OUTDIR}/${SAMPLE}.viral.summarized.csv"
echo "Bacteria/archaea list: ${OUTDIR}/${SAMPLE}.prok.summarized.csv"
echo "Readable version:       ${OUTDIR}/${SAMPLE}.viral.human.txt"
ls -lh "${OUTDIR}/"
