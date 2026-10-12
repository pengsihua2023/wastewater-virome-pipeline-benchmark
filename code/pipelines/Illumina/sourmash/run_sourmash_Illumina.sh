#!/bin/bash

#SBATCH --job-name=Sourmash_708_1-Illumina
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=256G
#SBATCH --time=5-04:00:00
#SBATCH --output=Sourmash_708-1-Illumina_S_%j.out
#SBATCH --error=Sourmash_708-1-Illumina_S_%j.err

# ============================================================================
# sourmash wastewater virome classification (single sample)
#   Run the NCBI Viruses viral database alone, scaled=50 at full resolution
#   classify phages and eukaryotic viruses together, no bacterial background, no filtering
#   Input: quality-controlled, host-removed paired-end reads
# Edited:26,29,30,31,46
# ============================================================================

set -euo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1

# ----------------------------------------------------------------------------
# 0. parameters -- change to your actual paths/samples
# ----------------------------------------------------------------------------
SAMPLE="sample_708_1"                          # sample name (output file prefix)

# Input reads (quality-controlled, host removed)
READS_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/illumina-708-genomes"
R1="${READS_DIR}/${SAMPLE}_01.fq.gz"
R2="${READS_DIR}/${SAMPLE}_02.fq.gz"

# Database
DB_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/sourmash_db"
VIRAL_DB="${DB_DIR}/ncbi-viruses-2025.01.k31.zip"
VIRAL_LIN="${DB_DIR}/ncbi-viruses.2025.01.lineages.csv"
TAX_SQLDB="${DB_DIR}/ncbi-viruses.2025.01.taxonomy.sqldb"   # see step 2; it is built automatically if missing

# sourmash parameters (must match the database: k=31, scaled=50)
KSIZE=31
SCALED=50
THRESHOLD_BP=100        # to catch very small viruses more sensitively, lower it to 50

# Output directory (placed under scratch, make sure it is within the container mount range and writable)
# if you want to change it, make sure the path is also inside /scratch, otherwise the container cannot write to it
OUTBASE="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/sourmash/Illumina-708/results_Illumina_708_1"
OUTDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "${OUTDIR}"

# Threads (aligned with the SLURM request)
THREADS=${SLURM_CPUS_PER_TASK:-2}

# ----------------------------------------------------------------------------
# 1. Load sourmash (Apptainer container: sourmash 4.9.4)
# ----------------------------------------------------------------------------
# .sif image location (you have pulled it into the containers directory)
SIF="/scratch/$USER/Meta-genome-data-analysis/Apptainer/containers/sourmash_4.9.4.sif"

# Wrap sourmash in a function: every later "sourmash xxx" call in the script then runs automatically inside the container
# --bind /scratch mounts the whole scratch into the container, so the container can read reads/databases and write output
# (reads, database and output directories are all under /scratch, so one mount covers everything)
sourmash() { apptainer exec --bind /scratch "${SIF}" sourmash "$@"; }

# If apptainer on your node needs module load first, uncomment the line below:
# module load apptainer

echo "==== $(date) | sample: ${SAMPLE} ===="
sourmash --version

# ----------------------------------------------------------------------------
# 2. Prepare the taxonomy sqldb (needed only once; skipped if it already exists)
# ----------------------------------------------------------------------------
if [[ ! -s "${TAX_SQLDB}" ]]; then
    echo "---- Building taxonomy sqldb ----"
    sourmash tax prepare \
        --taxonomy-csv "${VIRAL_LIN}" \
        -o "${TAX_SQLDB}" \
        -F sql
else
    echo "---- taxonomy sqldb already exists, skipping ----"
fi

# ----------------------------------------------------------------------------
# 3. SKETCH -- reads to signature (k=31, scaled=50, abund is mandatory)
# ----------------------------------------------------------------------------
SIG="${OUTDIR}/${SAMPLE}.sig.zip"
echo "---- sketch ----"
sourmash sketch dna \
    -p k=${KSIZE},scaled=${SCALED},abund \
    --name "${SAMPLE}" \
    -o "${SIG}" \
    "${R1}" "${R2}"

# check the signature parameters (to diagnose a k/scaled mismatch)
sourmash sig fileinfo "${SIG}"

# ----------------------------------------------------------------------------
# 4. GATHER -- against the viral database (single database, scaled=50 full resolution)
# ----------------------------------------------------------------------------
GATHER_CSV="${OUTDIR}/${SAMPLE}.viral.gather.csv"
UNASSIGNED="${OUTDIR}/${SAMPLE}.viral.unassigned.sig.zip"
echo "---- gather ----"

# gather returns non-zero when there is no viral match at all; this is normal for empty samples, so temporarily turn off -e
set +e
sourmash gather "${SIG}" "${VIRAL_DB}" \
    -k ${KSIZE} \
    --scaled ${SCALED} \
    --threshold-bp ${THRESHOLD_BP} \
    --estimate-ani-ci \
    -o "${GATHER_CSV}" \
    --output-unassigned "${UNASSIGNED}"
GATHER_RC=$?
set -e

if [[ ${GATHER_RC} -ne 0 && ! -s "${GATHER_CSV}" ]]; then
    echo "WARNING: sample ${SAMPLE}: no virus detected (gather found no match), skipping taxonomy"
    echo "==== $(date) | done (no viral hits) ===="
    exit 0
fi

# ----------------------------------------------------------------------------
# 5. TAX METAGENOME -- virus/phage classification profile (with full lineage, abundance-weighted)
# ----------------------------------------------------------------------------
echo "---- tax metagenome ----"
sourmash tax metagenome \
    --gather-csv "${GATHER_CSV}" \
    --taxonomy "${TAX_SQLDB}" \
    --output-base "${OUTDIR}/${SAMPLE}.viral" \
    --output-format csv_summary krona kreport human \
    --use-abundances \
    --rank species

# ----------------------------------------------------------------------------
# 6. Bacteria + archaea (GTDB rs226) -- independent gather
#    The viral database is scaled=50 and the GTDB k31 library is usually scaled=1000; mixing them in one gather would force both to the same value,
#    1000 with a loss of viral sensitivity, so GTDB is gathered separately for the same SIG.
#    scaled=50 sketch can be downsampled to 1000 (FracMinHash nesting), so it is compatible.
#    Short-read reads carry abundance information, so use --use-abundances.
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
        echo "---- tax metagenome (GTDB) ----"
        sourmash tax metagenome \
            --gather-csv "${GTDB_GATHER_CSV}" \
            --taxonomy "${GTDB_TAX_SQLDB}" \
            --output-base "${OUTDIR}/${SAMPLE}.prok" \
            --output-format csv_summary krona kreport human \
            --use-abundances \
            --rank species
        echo "Bacteria/archaea results: ${OUTDIR}/${SAMPLE}.prok.summarized.csv"
    else
        echo "WARNING: GTDB gather no match (${SAMPLE})"
    fi
fi

# ----------------------------------------------------------------------------
# Done
# ----------------------------------------------------------------------------
echo "==== $(date) | Done ===="
echo "Main result (viruses): ${OUTDIR}/${SAMPLE}.viral.summarized.csv"
echo "Bacteria/archaea:    ${OUTDIR}/${SAMPLE}.prok.summarized.csv"
echo "Readable version: ${OUTDIR}/${SAMPLE}.viral.human.txt"
ls -lh "${OUTDIR}/"
