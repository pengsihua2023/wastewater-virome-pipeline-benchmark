#!/bin/bash
#SBATCH --job-name=GO-Short_708-2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_sam0_Short_708_2_%j.out
#SBATCH --error=Viral_sam0_Short_708_2_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# Record the start time
START_TIME=$(date)
START_TIME_EPOCH=$(date +%s)

# Set a trap that records the end time on exit (use double quotes so that the variable is expanded correctly)
trap "END_TIME=\$(date); END_TIME_EPOCH=\$(date +%s); DURATION=\$((END_TIME_EPOCH - START_TIME_EPOCH)); echo ''; echo \"End time: \$END_TIME\"; echo \"Duration: \$((DURATION / 3600))h \$((DURATION % 3600 / 60))m \$((DURATION % 60))s\"" EXIT

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow (Short-Read Mode)"
echo "Data edited: 54, 64,68"
"=========================================="
echo "Start time: $START_TIME"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

# User's conda environment path
USER_CONDA_ENV="/home/$USER/.conda/envs/gottcha2_env"

# Activate conda environment
source activate "$USER_CONDA_ENV"

# Run GOTTCHA2 classification analysis
echo "🔬 2. Running GOTTCHA2 taxonomic profiling..."
# ============================================================
# GOTTCHA2 database settings  (NEW database — same one used by
# the PacBio run, dir "GOTTCHA2", 198G, contains viral signatures)
#   gottcha_db.species.fna.mmi
#   gottcha_db.species.fna.stats
#   gottcha_db.species.fna.tax.tsv
# -d -> .mmi index ; -ti -> .tax.tsv taxonomy
# ============================================================
DB_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/GOTTCHA2"
# IMPORTANT: pass the PREFIX (no .mmi). gottcha2 appends ".mmi" itself.
DB_PREFIX="${DB_DIR}/gottcha_db.species.fna"

# Sanity checks for DB and inputs (check the real files on disk)
for f in "${DB_PREFIX}.mmi" "${DB_PREFIX}.tax.tsv" sample_708_2_01.fq sample_708_2_02.fq; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required file not found or empty: $f" >&2
        exit 1
    fi
done

# Short-read (Illumina) paired-end profiling against the new DB
#   -xm sr : short-read minimap2 preset (NOT map-hifi)
gottcha2 profile \
  -i sample_708_2_01.fq sample_708_2_02.fq \
  -d "$DB_PREFIX" \
  -t 32 \
  -xm sr \
  -o result_shortread_sample_2_708
