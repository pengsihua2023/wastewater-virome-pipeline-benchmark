#!/bin/bash
#SBATCH --job-name=TaxTri_mina_708_2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=TaxTriage_Illumina_708_2_%j.out
#SBATCH --error=TaxTriage_Illumina_708_2_%j.err

set -euo pipefail

echo "=========================================="
echo "Job started at:"
date

echo "Running on node:"
hostname

echo "Current directory before cd:"
pwd
echo "=========================================="

# ------------------------------------------------------------------
# Java / conda environment (Nextflow needs Java 11-23; system default is Java 8)
# Activate the conda env that ships Java 17 + Nextflow.
# ------------------------------------------------------------------
module load Miniforge3/24.11.3-0 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
CONDA_ENV_PATH=/home/$USER/.conda/envs/nextflow_env
if [ -d "$CONDA_ENV_PATH" ]; then
    conda activate "$CONDA_ENV_PATH" || export PATH="$CONDA_ENV_PATH/bin:$PATH"
else
    conda activate nextflow_env
fi
# Login env may set PYTHONPATH/PYTHONHOME (user-py39) which breaks conda python3
unset PYTHONPATH PYTHONHOME

# Verify Java is 11-23 before continuing
JAVA_MAJOR=$(java -version 2>&1 | head -1 | grep -oE '"[0-9]+' | grep -oE '[0-9]+' | head -1)
echo "Java major version: ${JAVA_MAJOR:-unknown}"
if [ -z "$JAVA_MAJOR" ] || [ "$JAVA_MAJOR" -lt 11 ] || [ "$JAVA_MAJOR" -gt 23 ]; then
    echo "ERROR: Nextflow needs Java 11-23 but found '$JAVA_MAJOR'."
    echo "       'java' resolves to: $(which java)"
    echo "       Install Java in the env:  conda install -n nextflow_env openjdk=17"
    exit 1
fi

# TaxTriage pipeline directory
TAXTRIAGE_DIR=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/taxtriage

# Input, output, database and work directories; edited: 53,54,56
SAMPLESHEET=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/Illumina-708/samplesheet_Illumina_708.csv
OUTDIR=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/Illumina-708/results_illumina_708_2
KRAKEN2_DB=/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/k2_standard_20240605
NXF_WORK=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/work_illumina_708_2

# Nextflow / Apptainer/Singularity cache directory
export NXF_VER=26.04.3
export NXF_SYNTAX_PARSER=v1

export NXF_APPTAINER_CACHEDIR=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/apptainer_cache
export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/Meta-genome-data-analysis/Apptainer/TaxTriage/apptainer_cache

mkdir -p "$NXF_APPTAINER_CACHEDIR"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
mkdir -p "$NXF_WORK"
mkdir -p "$OUTDIR"

echo "=========================================="
echo "Checking paths..."

if [ ! -d "$TAXTRIAGE_DIR" ]; then
    echo "ERROR: TaxTriage directory not found:"
    echo "$TAXTRIAGE_DIR"
    exit 1
fi

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found:"
    echo "$SAMPLESHEET"
    exit 1
fi

if [ ! -d "$KRAKEN2_DB" ]; then
    echo "ERROR: Kraken2 database directory not found:"
    echo "$KRAKEN2_DB"
    exit 1
fi

echo "TaxTriage dir: $TAXTRIAGE_DIR"
echo "Samplesheet:   $SAMPLESHEET"
echo "Output dir:    $OUTDIR"
echo "Kraken2 DB:    $KRAKEN2_DB"
echo "Work dir:      $NXF_WORK"
echo "Cache dir:     $NXF_APPTAINER_CACHEDIR"
echo "=========================================="

# Enter the TaxTriage pipeline directory
cd "$TAXTRIAGE_DIR"

echo "Current directory after cd:"
pwd

echo "Nextflow version:"
nextflow -version

echo "Check container commands:"
which apptainer || true
which singularity || true

echo "=========================================="
echo "Running TaxTriage for Illumina data..."
echo "=========================================="

nextflow run main.nf \
  -profile singularity \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --db "$KRAKEN2_DB" \
  -work-dir "$NXF_WORK" \
#  -resume

echo "=========================================="
echo "Job finished at:"
date
echo "=========================================="