#!/bin/bash
#SBATCH --job-name=Pfiler_PacBio_708_2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_PacBio_708_2_%j.out
#SBATCH --error=Viral_PacBio_708_2_%j.err
#SBATCH --export=NONE

# === 1. Initialize Conda and activate nextflow_env (critical! contains Java 17 + apptainer) ===
# Initialize conda (try different common locations and methods) : edited line 272; remember to edit the config file too

# Method 1: Try to load conda module (common on HPC systems)
if command -v module &> /dev/null; then
    module load conda 2>/dev/null || module load anaconda 2>/dev/null || module load miniconda3 2>/dev/null || module load Miniforge3 2>/dev/null || true
fi

# Method 2: Try known Miniforge3 installation (from conda info output)
CONDA_INIT_SUCCESS=false
CONDA_BASE=""

# Check for Miniforge3 installation (common on HPC systems with EasyBuild)
if [ -f "/apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh" ]; then
    source "/apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh"
    CONDA_INIT_SUCCESS=true
    CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    echo "📍 Initialized conda from Miniforge3: $CONDA_BASE"
elif [ -f "/apps/eb/Miniforge3/24.11.3-0/bin/conda" ]; then
    eval "$(/apps/eb/Miniforge3/24.11.3-0/bin/conda shell.bash hook)"
    CONDA_INIT_SUCCESS=true
    CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    echo "📍 Initialized conda from Miniforge3: $CONDA_BASE"
# Try to find any Miniforge3 version
elif [ -d "/apps/eb/Miniforge3" ]; then
    # Find the latest or any Miniforge3 version
    LATEST_MINIFORGE=$(ls -td /apps/eb/Miniforge3/* 2>/dev/null | head -n 1)
    if [ -n "$LATEST_MINIFORGE" ] && [ -f "$LATEST_MINIFORGE/etc/profile.d/conda.sh" ]; then
        source "$LATEST_MINIFORGE/etc/profile.d/conda.sh"
        CONDA_INIT_SUCCESS=true
        CONDA_BASE="$LATEST_MINIFORGE"
        echo "📍 Initialized conda from Miniforge3: $CONDA_BASE"
    elif [ -n "$LATEST_MINIFORGE" ] && [ -f "$LATEST_MINIFORGE/bin/conda" ]; then
        eval "$($LATEST_MINIFORGE/bin/conda shell.bash hook)"
        CONDA_INIT_SUCCESS=true
        CONDA_BASE="$LATEST_MINIFORGE"
        echo "📍 Initialized conda from Miniforge3: $CONDA_BASE"
    fi
fi

# Method 3: Infer conda base from known environment path (if Method 2 failed)
if [ "$CONDA_INIT_SUCCESS" = false ] && [ -d "/scratch/$USER/conda/envs/nextflow_env" ]; then
    # Environment path: /scratch/$USER/conda/envs/nextflow_env
    # Try to find conda base that manages this environment
    echo "📍 Found nextflow_env at: /scratch/$USER/conda/envs/nextflow_env"
    
    # Try common conda base locations
    for potential_base in "/scratch/$USER/conda" "/apps/eb/Miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
        if [ -f "$potential_base/etc/profile.d/conda.sh" ]; then
            source "$potential_base/etc/profile.d/conda.sh"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$potential_base"
            echo "   ✅ Initialized using: $CONDA_BASE/etc/profile.d/conda.sh"
            break
        elif [ -f "$potential_base/bin/conda" ]; then
            eval "$($potential_base/bin/conda shell.bash hook)"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$potential_base"
            echo "   ✅ Initialized using: $CONDA_BASE/bin/conda shell.bash hook"
            break
        fi
    done
fi

# Method 4: Try to find conda in PATH (if already initialized)
if [ "$CONDA_INIT_SUCCESS" = false ] && command -v conda &> /dev/null; then
    CONDA_BASE_FOUND=$(conda info --base 2>/dev/null)
    if [ -n "$CONDA_BASE_FOUND" ]; then
        if [ -f "$CONDA_BASE_FOUND/etc/profile.d/conda.sh" ]; then
            source "$CONDA_BASE_FOUND/etc/profile.d/conda.sh"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$CONDA_BASE_FOUND"
        elif [ -f "$CONDA_BASE_FOUND/bin/conda" ]; then
            eval "$($CONDA_BASE_FOUND/bin/conda shell.bash hook)"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$CONDA_BASE_FOUND"
        fi
    fi
fi

# Method 5: Try standard locations (if previous methods failed)
if [ "$CONDA_INIT_SUCCESS" = false ]; then
    for conda_base in "$HOME/anaconda3" "$HOME/miniconda3" "/opt/anaconda3" "/opt/miniconda3"; do
        if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
            source "$conda_base/etc/profile.d/conda.sh"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$conda_base"
            break
        elif [ -f "$conda_base/bin/conda" ]; then
            eval "$($conda_base/bin/conda shell.bash hook)"
            CONDA_INIT_SUCCESS=true
            CONDA_BASE="$conda_base"
            break
        fi
    done
fi

# Final check
if [ "$CONDA_INIT_SUCCESS" = false ] || ! command -v conda &> /dev/null; then
    echo "❌ Error: Cannot find or initialize conda."
    echo ""
    echo "   Debugging information:"
    echo "   - Environment path exists: $([ -d /scratch/$USER/conda/envs/nextflow_env ] && echo 'yes' || echo 'no')"
    echo "   - Miniforge3 exists: $([ -d /apps/eb/Miniforge3/24.11.3-0 ] && echo 'yes' || echo 'no')"
    echo "   - Miniforge3 init script: $([ -f /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh ] && echo 'yes' || echo 'no')"
    echo "   - Miniforge3 binary: $([ -f /apps/eb/Miniforge3/24.11.3-0/bin/conda ] && echo 'yes' || echo 'no')"
    echo "   - Conda base (/scratch/$USER/conda) exists: $([ -d /scratch/$USER/conda ] && echo 'yes' || echo 'no')"
    echo ""
    echo "   Solution options:"
    echo "   1. Try loading Miniforge3 module: module load Miniforge3"
    echo "   2. Check if conda is in PATH: which conda"
    echo "   3. Find conda base: conda info --base"
    echo "   4. Manually source: source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh"
    exit 1
fi

# Activate nextflow_env environment
# Set CONDA_ENVS_PATH to include custom environment directory if environment is not in default location
# Use the HOME conda env explicitly. The same-named env under /scratch has a
# BROKEN python (cannot load its own stdlib -> "init_fs_encoding / can't find
# encoding"), which broke the abundance step. The HOME env has a working python
# with pandas/numpy. NOTE: nextflow & java for the main analysis come from
# `module load` (below), NOT from this conda env, so switching envs is safe.
ENV_PATH="/home/$USER/.conda/envs/nextflow_env"
ENV_NAME="nextflow_env"

if [ -d "$ENV_PATH" ]; then
    # Environment exists in custom location
    # Try to add the parent directory to conda envs path
    export CONDA_ENVS_PATH="/scratch/$USER/conda/envs:$CONDA_ENVS_PATH"
    
    # Try to create symlink in conda base/envs directory for easier activation
    # Find conda base directory
    if [ -n "$CONDA_BASE" ]; then
        CONDA_BASE_DIR="$CONDA_BASE"
    elif command -v conda &> /dev/null; then
        CONDA_BASE_DIR=$(conda info --base 2>/dev/null)
    elif [ -d "/apps/eb/Miniforge3/24.11.3-0" ]; then
        CONDA_BASE_DIR="/apps/eb/Miniforge3/24.11.3-0"
    else
        CONDA_BASE_DIR=""
    fi
    
    # Create symlink if conda base found and link doesn't exist
    if [ -n "$CONDA_BASE_DIR" ] && [ -d "$CONDA_BASE_DIR/envs" ]; then
        SYMLINK_PATH="$CONDA_BASE_DIR/envs/$ENV_NAME"
        if [ ! -e "$SYMLINK_PATH" ]; then
            echo "🔗 Creating symlink: $SYMLINK_PATH -> $ENV_PATH"
            ln -s "$ENV_PATH" "$SYMLINK_PATH" 2>/dev/null && echo "   ✅ Symlink created" || echo "   ⚠️  Failed to create symlink (may need permissions)"
        elif [ -L "$SYMLINK_PATH" ]; then
            CURRENT_TARGET=$(readlink "$SYMLINK_PATH")
            if [ "$CURRENT_TARGET" != "$ENV_PATH" ]; then
                echo "🔗 Updating symlink: $SYMLINK_PATH -> $ENV_PATH"
                rm "$SYMLINK_PATH" 2>/dev/null
                ln -s "$ENV_PATH" "$SYMLINK_PATH" 2>/dev/null && echo "   ✅ Symlink updated" || echo "   ⚠️  Failed to update symlink"
            fi
        fi
    fi
    
    # Activate by FULL PATH (not by name) so conda cannot resolve "nextflow_env"
    # to the broken /scratch env of the same name.
    if ! conda activate "$ENV_PATH" 2>/dev/null; then
        # If normal activation fails, try with full path
        echo "📍 Environment found at: $ENV_PATH"
        echo "   Attempting activation with full path..."
        
        # Try direct activation using conda env list or source activation script
        if [ -f "$ENV_PATH/bin/activate" ]; then
            source "$ENV_PATH/bin/activate"
            echo "   ✅ Activated using full path: $ENV_PATH"
        elif [ -f "$ENV_PATH/etc/conda/activate.d/activate.sh" ]; then
            source "$ENV_PATH/etc/conda/activate.d/activate.sh"
            echo "   ✅ Activated using activation script"
        else
            # Last resort: set PATH manually
            export PATH="$ENV_PATH/bin:$PATH"
            export CONDA_DEFAULT_ENV="nextflow_env"
            export CONDA_PREFIX="$ENV_PATH"
            echo "   ✅ Activated by setting PATH and environment variables"
        fi
        
        # Verify activation
        if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "nextflow_env" ]; then
            echo "⚠️  Warning: Environment may not be fully activated. Continuing anyway..."
        fi
    else
        echo "✅ Activated nextflow_env successfully"
    fi
else
    # Try activation by full path if the expected dir wasn't found above
    conda activate "$ENV_PATH" || {
        echo "❌ Error: Cannot activate nextflow_env. Environment not found."
        echo ""
        echo "   Debugging information:"
        echo "   - Expected path: $ENV_PATH"
        echo "   - Path exists: $([ -d "$ENV_PATH" ] && echo 'yes' || echo 'no')"
        echo "   - Available environments:"
        conda env list 2>/dev/null || echo "   Cannot list environments"
        echo ""
        echo "   Solution: Please ensure the environment exists at $ENV_PATH"
        echo "   Or create it with: conda create -n nextflow_env python=3.9 nextflow java=17 -p $ENV_PATH"
        exit 1
    }
fi

# === 2. Verify Java version (must be 11-23) ===
JAVA_VERSION=$(java -version 2>&1 | head -n 1 | sed -E 's/.*version "([0-9]+)\..*/\1/')
if [ -z "$JAVA_VERSION" ] || [ "$JAVA_VERSION" -lt 11 ] || [ "$JAVA_VERSION" -gt 23 ]; then
    echo "⚠️  Warning: Java version may be incompatible (detected: $JAVA_VERSION, required: 11-23)"
    echo "   Attempting to load Java module or use conda Java..."
    # Try to load Java module if available
    module load java/17 2>/dev/null || module load java/11 2>/dev/null || true
fi

# === 3. Load Nextflow (if not in environment) ===
module load nextflow 2>/dev/null || true

# === 4. Set Apptainer cache ===
export APPTAINER_CACHEDIR=$HOME/.apptainer_cache
export NXF_SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR

# === 5. Verify environment (optional, for debugging) ===
echo "=== Environment Verification ==="
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "Java:"
java -version
echo "Apptainer:"
apptainer --version
echo "Nextflow:"
nextflow -v

# === 6. Run nf-core/taxprofiler workflow (PacBio HiFi long-read data) ===
echo ""
echo "========================================"
echo "▶️  Starting TaxProfiler Analysis"
echo "========================================"

nextflow run nf-core/taxprofiler \
  -r 1.2.0 \
  -profile apptainer \
  -c nextflow_pacbio.config \
  -resume

# Check if workflow completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✅ TaxProfiler Analysis Completed"
    echo "========================================"
    
    # === 7. Automatic viral abundance calculation (RPM & RPKM) ===
    echo ""
    echo "========================================"
    echo "🧮 Auto-calculating Viral Abundance (RPM & RPKM)"
    echo "========================================"
    
    # Get output directory (read from config file)
    OUTDIR=$(grep "outdir" nextflow_pacbio.config | head -1 | sed "s/.*= *['\"]//;s/['\"].*//")
    
    if [ -z "$OUTDIR" ]; then
        OUTDIR="results_pacbio_sample_2"
    fi
    
    echo "Output directory: $OUTDIR"
    
    # Run abundance calculation script (long-read specific)
    # Note: Long-read data doesn't use Bracken, calculates abundance directly from Kraken2 report
    if [ -f "batch_calculate_abundance_longread_en.sh" ]; then
        # The login env sets PYTHONPATH/PYTHONHOME (user-py39) which breaks the
        # conda Python used by the abundance step. Clear them before calling it.
        unset PYTHONPATH PYTHONHOME

        # Export the conda environment variables so that child scripts can detect the conda environment
        export CONDA_PREFIX
        export CONDA_DEFAULT_ENV
        export CONDA_ENVS_PATH

        # If the environment is already activated, explicitly set the Python path.
        # Fallback points at the HOME env (working python), NOT the /scratch env
        # (broken python). The abundance script also re-validates the interpreter
        # and auto-skips any candidate whose pandas/numpy import fails.
        if [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/bin/python" ]; then
            export PYTHON_CMD="$CONDA_PREFIX/bin/python"
            echo "   Setting PYTHON_CMD: $PYTHON_CMD"
        elif [ -f "$ENV_PATH/bin/python" ]; then
            export PYTHON_CMD="$ENV_PATH/bin/python"
            export CONDA_PREFIX="$ENV_PATH"
            export CONDA_DEFAULT_ENV="nextflow_env"
            echo "   Setting PYTHON_CMD: $PYTHON_CMD"
        fi
        
        bash batch_calculate_abundance_longread_en.sh "$OUTDIR"
        
        if [ $? -eq 0 ]; then
            echo ""
            echo "========================================"
            echo "🎉 All Tasks Completed!"
            echo "========================================"
            echo "📊 Classification results: $OUTDIR/kraken2/"
            echo "📈 Abundance results: $OUTDIR/abundance/"
            echo "📋 Comprehensive report: $OUTDIR/multiqc/multiqc_report.html"
            echo "========================================"
        else
            echo "⚠️  Warning: Abundance calculation failed, but main analysis completed"
        fi
    else
        echo "⚠️  Warning: batch_calculate_abundance_longread_en.sh script not found"
        echo "    You can run manually: bash batch_calculate_abundance_longread_en.sh $OUTDIR"
    fi
else
    echo ""
    echo "========================================"
    echo "❌ TaxProfiler Analysis Failed"
    echo "========================================"
    echo "Please check log files for detailed information"
    exit 1
fi

