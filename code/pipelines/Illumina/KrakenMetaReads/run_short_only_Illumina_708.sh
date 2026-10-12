#!/bin/bash
#SBATCH --job-name=Kraken-Illumina_708_0
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=15-00:00:00
#SBATCH --output=Short_708_0_%j.out
#SBATCH --error=Short_708_0_%j.err

set -o pipefail

cd "$SLURM_SUBMIT_DIR" || exit 1

# Edited:16
OUTDIR_SHORT="results_short_Illumina_708_0"
KRAKEN2_DB="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/k2_standard_20240605"
PROJECT_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/KrakenMetaReads-nf"
CONDA_ENV="/scratch/$USER/conda/envs/nextflow_env"
NXF_CONDA_CACHE="/scratch/$USER/nextflow_conda_cache"

echo "=========================================="
echo "🧬  Short-Read Metagenome Workflow"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node: ${SLURM_NODELIST:-N/A}"
echo "Submit dir: $SLURM_SUBMIT_DIR"
echo "Output dir: $OUTDIR_SHORT"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Load Java module or use Java from conda environment
echo "☕ Loading Java module..."
JAVA_LOADED=false

# Try common Java module names first.
for java_module in Java/21 Java/17 Java/11 Java/21.0.2 Java/17.0.8 Java/11.0.20; do
    if module load "$java_module" 2>/dev/null; then
        echo "✅ Loaded Java module: $java_module"
        JAVA_LOADED=true
        break
    fi
done

# Fall back to Java inside the conda environment.
if [ "$JAVA_LOADED" = false ]; then
    echo "⚠️  Java module not loaded. Checking conda environment..."
    if [ -x "$CONDA_ENV/bin/java" ]; then
        export JAVA_HOME="$CONDA_ENV"
        export PATH="$JAVA_HOME/bin:$PATH"
        echo "✅ Using Java from conda environment: $CONDA_ENV/bin/java"
    else
        echo "❌ Java 11-23 not found. Please install Java in the Nextflow conda environment:"
        echo "   conda activate $CONDA_ENV"
        echo "   conda install -c conda-forge openjdk=17"
        exit 1
    fi
fi

# Verify Java version.
JAVA_VERSION=$(java -version 2>&1 | head -n 1)
JAVA_VERSION_STR=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
echo "   Java: $JAVA_VERSION"

if [ -z "$JAVA_VERSION_STR" ]; then
    echo "❌ Could not determine Java version."
    exit 1
fi

if echo "$JAVA_VERSION_STR" | grep -q "^1\."; then
    JAVA_MAJOR=$(echo "$JAVA_VERSION_STR" | awk -F '.' '{print $2}')
else
    JAVA_MAJOR=$(echo "$JAVA_VERSION_STR" | awk -F '.' '{print $1}')
fi

if [ "$JAVA_MAJOR" -lt 11 ] || [ "$JAVA_MAJOR" -gt 23 ]; then
    echo "❌ Java version $JAVA_VERSION_STR is not supported. Nextflow requires Java 11-23."
    exit 1
fi

echo "✅ Java version $JAVA_MAJOR is compatible with Nextflow"

# Verify tools
echo ""
echo "🧪 2. Verifying tools..."
if command -v nextflow &>/dev/null; then
    echo "✅ Nextflow: $(which nextflow)"
else
    echo "❌ Nextflow not found in PATH after activating $CONDA_ENV"
    exit 1
fi

if command -v apptainer &>/dev/null; then
    echo "✅ Apptainer: $(which apptainer)"
elif command -v singularity &>/dev/null; then
    echo "✅ Singularity: $(which singularity)"
else
    echo "❌ Apptainer/Singularity not found"
    exit 1
fi

echo ""
echo "ℹ️  Workflow: Short-read only (MEGAHIT + SPAdes + Bowtie2 + Kraken2)"
echo ""

# Set short-read samplesheet path.
if [ -f "samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="samplesheet_short.csv"
elif [ -f "$PROJECT_DIR/samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="$PROJECT_DIR/samplesheet_short.csv"
else
    echo "❌ Short-read samplesheet not found!"
    echo "   Checked: ./samplesheet_short.csv"
    echo "   Checked: $PROJECT_DIR/samplesheet_short.csv"
    exit 1
fi

# Verify databases
echo "🗄️ 3. Verifying databases..."
if [ -d "$KRAKEN2_DB" ]; then
    echo "✅ Kraken2 database: $KRAKEN2_DB"
else
    echo "❌ Kraken2 database not found: $KRAKEN2_DB"
    exit 1
fi
echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "$SAMPLESHEET_SHORT" ]; then
    echo "✅ Short-read samplesheet: $SAMPLESHEET_SHORT"
    SHORT_SAMPLES=$(tail -n +2 "$SAMPLESHEET_SHORT" | sed '/^[[:space:]]*$/d' | wc -l)
    echo "   📊 Found $SHORT_SAMPLES short-read samples"
else
    echo "❌ Samplesheet not found: $SAMPLESHEET_SHORT"
    exit 1
fi
echo ""

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ -d "$OUTDIR_SHORT" ]; then
    echo "Removing previous short-read results: $OUTDIR_SHORT"
    rm -rf "$OUTDIR_SHORT"
fi
echo ""

# Set Apptainer/Singularity bind paths.
# Binding /scratch is safer because the samplesheet, reads, database, and work directories are all likely under /scratch.
export APPTAINER_BIND="/scratch:/scratch"
export SINGULARITY_BIND="/scratch:/scratch"

# Set Nextflow/conda environment variables.
echo "🔧 Setting up conda environment for Nextflow processes..."
export CONDA_BASE
export PATH="$CONDA_BASE/bin:$PATH"

if [ -x "$CONDA_BASE/bin/conda" ]; then
    CONDA_PATH="$CONDA_BASE/bin/conda"
elif command -v conda &>/dev/null && ! type conda | grep -q "function"; then
    CONDA_PATH=$(which conda)
else
    CONDA_PATH="$CONDA_BASE/bin/conda"
fi

mkdir -p "$NXF_CONDA_CACHE"
export NXF_CONDA_CACHEDIR="$NXF_CONDA_CACHE"
export NXF_CONDA_PATH="$CONDA_PATH"

echo "✅ Conda base: $CONDA_BASE"
echo "✅ Conda path: $CONDA_PATH"
echo "✅ Nextflow conda cache: $NXF_CONDA_CACHEDIR"
echo "✅ APPTAINER_BIND: $APPTAINER_BIND"
echo ""

# Run workflow
echo "🚀 6. Running short-read workflow..."
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control"
echo "   2. MEGAHIT and metaSPAdes assembly"
echo "   3. Bowtie2 mapping and abundance"
echo "   4. Kraken2 classification"
echo "   5. Comprehensive report generation"
echo ""

nextflow run metagenome_hybrid_workflow.nf \
    -c metagenome_hybrid_workflow.config \
    --input_short "$SAMPLESHEET_SHORT" \
    --outdir_short "$OUTDIR_SHORT" \
    --kraken2_db "$KRAKEN2_DB"

WORKFLOW_EXIT_CODE=$?

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ "$WORKFLOW_EXIT_CODE" -eq 0 ]; then
    echo "✅ Workflow completed successfully!"

    if [ -d "$OUTDIR_SHORT" ]; then
        echo ""
        echo "📊 Short-read results ($OUTDIR_SHORT/):"

        if [ -d "$OUTDIR_SHORT/fastp" ]; then
            FASTP_COUNT=$(find "$OUTDIR_SHORT/fastp" -name "*.html" 2>/dev/null | wc -l)
            echo "  ✅ fastp reports: $FASTP_COUNT files"
        fi

        if [ -d "$OUTDIR_SHORT/abundance_megahit" ]; then
            ABUND_M=$(find "$OUTDIR_SHORT/abundance_megahit" -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ MEGAHIT abundance: $ABUND_M files"
        fi

        if [ -d "$OUTDIR_SHORT/abundance_spades" ]; then
            ABUND_S=$(find "$OUTDIR_SHORT/abundance_spades" -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ SPAdes abundance: $ABUND_S files"
        fi

        if [ -d "$OUTDIR_SHORT/kraken2_megahit" ]; then
            echo "  ✅ Kraken2 MEGAHIT classification"
        fi

        if [ -d "$OUTDIR_SHORT/kraken2_spades" ]; then
            echo "  ✅ Kraken2 SPAdes classification"
        fi

        if [ -d "$OUTDIR_SHORT/merged_reports" ]; then
            MERGED=$(find "$OUTDIR_SHORT/merged_reports" -name "*_merged_report.csv" 2>/dev/null | wc -l)
            CONSENSUS=$(find "$OUTDIR_SHORT/merged_reports" -name "*_virus_consensus.txt" 2>/dev/null | wc -l)
            echo "  ✅ Merged reports: $MERGED files"
            if [ "$CONSENSUS" -gt 0 ]; then
                echo "  ✅ Virus consensus analysis: $CONSENSUS files ⭐"
            fi
        fi

        echo ""
        echo "📋 Summary:"
        echo "  Total files: $(find "$OUTDIR_SHORT" -type f 2>/dev/null | wc -l)"
        echo ""
        echo "📁 Results location: $OUTDIR_SHORT/"
    else
        echo "❌ Results directory not found: $OUTDIR_SHORT"
    fi
else
    echo "❌ Workflow failed with exit code: $WORKFLOW_EXIT_CODE"
    echo "🔍 Check the error log for details: Short_Read_sam0_${SLURM_JOB_ID}.err"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="

exit "$WORKFLOW_EXIT_CODE"
