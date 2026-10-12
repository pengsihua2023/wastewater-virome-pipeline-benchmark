#!/bin/bash
#SBATCH --job-name=PacBio_Kraks_Sam2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=96:00:00
#SBATCH --output=PacBio_KrakenMetaReads_sam2_%j.out
#SBATCH --error=PacBio_KrakenMetaReads_sam2_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# -----------------------------------------------------------------------
# OUTPUT DIRECTORY — change only this one variable to switch sample runs
# Edited:17
# -----------------------------------------------------------------------
OUTDIR="results_pacbio_sample-2"
# -----------------------------------------------------------------------

echo "=========================================="
echo "PacBio HiFi Metagenome Workflow"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Output dir: $OUTDIR"
echo ""

# Load conda environment
echo "1. Setting up environment..."
module load Miniforge3/24.11.3-0
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate /scratch/$USER/conda/envs/nextflow_env

# Load Java (Nextflow requires Java 11-23)
echo "Loading Java module..."
for java_version in 21 17 11; do
    if module load Java/${java_version} 2>/dev/null; then
        echo "Loaded Java ${java_version}"
        break
    fi
done

if ! java -version &>/dev/null; then
    if [ -f "/scratch/$USER/conda/envs/nextflow_env/bin/java" ]; then
        export JAVA_HOME="/scratch/$USER/conda/envs/nextflow_env"
        export PATH="$JAVA_HOME/bin:$PATH"
        echo "Using Java from conda environment"
    else
        echo "ERROR: Java 11+ not found."
        echo "  conda activate /scratch/$USER/conda/envs/nextflow_env"
        echo "  conda install -c conda-forge openjdk=11"
        exit 1
    fi
fi

# Validate Java version (Nextflow requires 11-23)
JAVA_VERSION_STR=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
if echo "$JAVA_VERSION_STR" | grep -q "^1\."; then
    echo "ERROR: Java 8 ($JAVA_VERSION_STR) is not supported. Nextflow requires Java 11-23."
    exit 1
else
    JAVA_MAJOR=$(echo "$JAVA_VERSION_STR" | awk -F '.' '{print $1}')
    if [ "$JAVA_MAJOR" -lt 11 ] || [ "$JAVA_MAJOR" -gt 23 ]; then
        echo "ERROR: Java $JAVA_MAJOR is not supported. Nextflow requires Java 11-23."
        exit 1
    fi
    echo "Java $JAVA_MAJOR is compatible with Nextflow"
fi

# Verify required tools
echo ""
echo "2. Verifying tools..."
NEXTFLOW_PATH=$(which nextflow)
echo "Nextflow: $NEXTFLOW_PATH"

if command -v apptainer &>/dev/null; then
    echo "Apptainer: $(which apptainer)"
elif command -v singularity &>/dev/null; then
    echo "Singularity: $(which singularity)"
else
    echo "ERROR: Apptainer/Singularity not found"
    exit 1
fi

if conda env list | grep -q "viralFlye_env"; then
    echo "viralFlye_env: Found"
else
    echo "ERROR: viralFlye_env not found!"
    echo "  conda create -n viralFlye_env python=3.8"
    echo "  conda activate viralFlye_env"
    echo "  git clone https://github.com/Dmitry-Antipov/viralFlye.git"
    echo "  cd viralFlye && pip install -r requirements.txt"
    exit 1
fi

echo ""
echo "Workflow execution environment:"
echo "  - Read type:     PacBio HiFi (--pacbio-hifi)"
echo "  - Minimap2 preset: asm20 (HiFi-optimized)"
echo "  - Assembly (metaFlye):    Apptainer"
echo "  - Mapping (Minimap2, Samtools): Conda"
echo "  - Viral identification (viralFlye): Conda"
echo "  - Classification (Kraken2): Conda"
echo ""

# Database and samplesheet paths
KRAKEN2_DB="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

if [ -f "samplesheet_pacbio.csv" ]; then
    SAMPLESHEET_PACBIO="samplesheet_pacbio.csv"
else
    echo "ERROR: PacBio samplesheet not found: samplesheet_pacbio.csv"
    exit 1
fi

# Verify databases
echo "3. Verifying databases..."
if [ -d "$KRAKEN2_DB" ]; then
    echo "Kraken2 database: $KRAKEN2_DB"
else
    echo "ERROR: Kraken2 database not found: $KRAKEN2_DB"
    exit 1
fi
echo ""

# Verify input files
echo "4. Verifying input files..."
PACBIO_SAMPLES=$(tail -n +2 "$SAMPLESHEET_PACBIO" | grep -c "." || true)
echo "PacBio samplesheet: $SAMPLESHEET_PACBIO ($PACBIO_SAMPLES samples)"

# Verify each data file listed in the samplesheet exists
while IFS=',' read -r sample fastq_long; do
    [ "$sample" = "sample" ] && continue  # skip header
    [ -z "$fastq_long" ] && continue
    if [ ! -f "$fastq_long" ]; then
        echo "ERROR: Input file not found: $fastq_long"
        exit 1
    fi
    echo "  OK: $sample -> $fastq_long"
done < "$SAMPLESHEET_PACBIO"
echo ""

# Clean previous results
echo "5. Cleaning previous results..."
if [ -d "$OUTDIR" ]; then
    echo "Removing previous PacBio results: $OUTDIR"
    rm -rf "$OUTDIR"
fi
echo ""

# Singularity bind paths
export SINGULARITY_BIND="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases:/databases"

# Export conda paths for Nextflow sub-processes
export CONDA_BASE=$(conda info --base)
export PATH="$CONDA_BASE/bin:$PATH"
CONDA_PATH="$CONDA_BASE/bin/conda"
if [ ! -f "$CONDA_PATH" ]; then
    echo "WARNING: conda binary not found at $CONDA_PATH"
fi
export NXF_CONDA_CACHEDIR="$CONDA_BASE"
export NXF_CONDA_PATH="$CONDA_PATH"
echo "Conda base: $CONDA_BASE"
echo "Conda path: $CONDA_PATH"
echo ""

# Run workflow
echo "6. Running PacBio HiFi long-read workflow..."
echo ""
echo "Pipeline steps:"
echo "  1. metaFlye assembly (--pacbio-hifi, metagenome mode)"
echo "  2. viralFlye viral contig identification"
echo "     - Linear viral contigs"
echo "     - Circular viral contigs"
echo "  3. Minimap2 mapping and abundance (RPM/RPKM)"
echo "     - metaFlye all contigs          [asm20 preset]"
echo "     - viralFlye linear viral contigs [asm20 preset]"
echo "     - viralFlye circular viral contigs [asm20 preset]"
echo "  4. Kraken2 classification for all contig sets"
echo ""

nextflow run metagenome_hybrid_workflow.nf \
    -c metagenome_pacbio_workflow.config \
    --input_long "$SAMPLESHEET_PACBIO" \
    --outdir_long "$OUTDIR" \
    --kraken2_db "$KRAKEN2_DB" \
    --long_read_type pacbio-hifi

WORKFLOW_EXIT_CODE=$?

# Report results
echo ""
echo "=========================================="
echo "Workflow Results"
echo "=========================================="

if [ $WORKFLOW_EXIT_CODE -eq 0 ]; then
    echo "Workflow completed successfully!"

    if [ -d "$OUTDIR" ]; then
        echo ""
        echo "PacBio results ($OUTDIR/):"

        for dir_label in \
            "abundance_flye:metaFlye abundance" \
            "viralflye:viralFlye viral contigs" \
            "abundance_viralflye_linear:viralFlye linear viral abundance" \
            "abundance_viralflye_circular:viralFlye circular viral abundance" \
            "kraken2_flye:Kraken2 metaFlye classification" \
            "kraken2_viralflye_linear:Kraken2 linear viral classification" \
            "kraken2_viralflye_circular:Kraken2 circular viral classification"
        do
            dir="${dir_label%%:*}"
            label="${dir_label##*:}"
            full_path="$OUTDIR/${dir}"
            if [ -d "$full_path" ]; then
                count=$(find "$full_path" -type f 2>/dev/null | wc -l)
                echo "  $label: $count files"
            fi
        done

        echo ""
        TOTAL=$(find "$OUTDIR" -type f 2>/dev/null | wc -l)
        echo "Total output files: $TOTAL"
        echo "Results location:   $OUTDIR/"
    else
        echo "WARNING: Results directory not found: $OUTDIR"
    fi

else
    echo "ERROR: Workflow failed with exit code $WORKFLOW_EXIT_CODE"
    echo "Check the error log: PacBio_KrakenMetaReads_sam2_${SLURM_JOB_ID}.err"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
