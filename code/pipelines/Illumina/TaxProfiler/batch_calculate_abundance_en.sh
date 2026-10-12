#!/bin/bash
# Batch calculate viral abundance for all samples (RPM & RPKM)

# === Detect and use correct Python interpreter ===
# Prefer the Python from the conda environment so that pandas and numpy are available

PYTHON_CMD=""

# Method 0: if the parent script has already set PYTHON_CMD, use it directly (highest priority)
if [ -n "$PYTHON_CMD" ] && [ -f "$PYTHON_CMD" ]; then
    echo "✅ Using Python from parent script: $PYTHON_CMD"
# Method 1: check CONDA_PREFIX (most reliable)
elif [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/bin/python" ]; then
    PYTHON_CMD="$CONDA_PREFIX/bin/python"
    echo "✅ Using conda Python from CONDA_PREFIX: $PYTHON_CMD"
# Method 2: check known custom environment paths
elif [ -f "/scratch/$USER/conda/envs/nextflow_env/bin/python" ]; then
    PYTHON_CMD="/scratch/$USER/conda/envs/nextflow_env/bin/python"
    echo "✅ Using conda Python from custom path: $PYTHON_CMD"
# Method 3: check symlinks under conda base/envs
elif [ -L "/apps/eb/Miniforge3/24.11.3-0/envs/nextflow_env" ] && [ -f "/apps/eb/Miniforge3/24.11.3-0/envs/nextflow_env/bin/python" ]; then
    PYTHON_CMD="/apps/eb/Miniforge3/24.11.3-0/envs/nextflow_env/bin/python"
    echo "✅ Using conda Python from symlink: $PYTHON_CMD"
# Method 4: try to find python in PATH (if the conda environment is activated)
elif [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" = "nextflow_env" ] && command -v python &> /dev/null; then
    PYTHON_PATH=$(which python)
    # Verify that this python is in the conda environment
    if [[ "$PYTHON_PATH" == *"conda"* ]] || [[ "$PYTHON_PATH" == *"nextflow_env"* ]]; then
        PYTHON_CMD="python"
        echo "✅ Using conda Python from PATH: $PYTHON_PATH"
    fi
fi

# If still not found, try other methods
if [ -z "$PYTHON_CMD" ]; then
    # Try to find all possible conda environment paths
    for env_path in \
        "/scratch/$USER/conda/envs/nextflow_env/bin/python" \
        "$HOME/.conda/envs/nextflow_env/bin/python" \
        "/opt/conda/envs/nextflow_env/bin/python"; do
        if [ -f "$env_path" ]; then
            PYTHON_CMD="$env_path"
            echo "✅ Found conda Python: $PYTHON_CMD"
            break
        fi
    done
fi

# Last fallback
if [ -z "$PYTHON_CMD" ]; then
    if command -v python &> /dev/null; then
        PYTHON_CMD="python"
        PYTHON_PATH=$(which python)
        echo "⚠️  Warning: Using Python from PATH: $PYTHON_PATH"
        echo "   This may not be the correct conda environment Python!"
    elif command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
        PYTHON_PATH=$(which python3)
        echo "⚠️  Warning: Using python3 from PATH: $PYTHON_PATH"
        echo "   This may not be the correct conda environment Python!"
    else
        echo "❌ Error: No Python interpreter found"
        exit 1
    fi
fi

# Verify that Python is available
if ! command -v "$PYTHON_CMD" &> /dev/null && [ ! -f "$PYTHON_CMD" ]; then
    echo "❌ Error: Python interpreter not found: $PYTHON_CMD"
    exit 1
fi

# Show Python information
echo "   Python version: $($PYTHON_CMD --version 2>&1)"
echo "   Python location: $($PYTHON_CMD -c 'import sys; print(sys.executable)' 2>/dev/null || echo "$PYTHON_CMD")"
echo ""

# Quickly check whether pandas and numpy are available
echo "🔍 Checking for pandas and numpy..."
if ! "$PYTHON_CMD" -c "import pandas, numpy" 2>/dev/null; then
    echo "❌ Error: pandas or numpy not available in Python environment"
    echo ""
    echo "   Debug information:"
    echo "   - Python: $($PYTHON_CMD --version 2>&1)"
    echo "   - Python path: $($PYTHON_CMD -c 'import sys; print(sys.executable)' 2>/dev/null || echo "$PYTHON_CMD")"
    echo "   - Python sys.path:"
    "$PYTHON_CMD" -c "import sys; [print(f'     {p}') for p in sys.path]" 2>/dev/null || echo "     (cannot get sys.path)"
    echo ""
    echo "   Solution: Install pandas and numpy in the conda environment:"
    echo "   conda activate nextflow_env"
    echo "   conda install pandas numpy -c conda-forge -y"
    echo ""
    echo "   Or if using custom path:"
    echo "   /scratch/$USER/conda/envs/nextflow_env/bin/pip install pandas numpy"
    exit 1
else
    echo "✅ pandas and numpy are available"
    "$PYTHON_CMD" -c "import pandas, numpy; print(f'   pandas: {pandas.__version__}, numpy: {numpy.__version__}')" 2>/dev/null
fi
echo ""

# Usage instructions
usage() {
    cat << EOF
Usage: bash batch_calculate_abundance_en.sh <results_dir> [output_dir]

Parameters:
    results_dir   TaxProfiler output results directory (e.g. results_viral_short)
    output_dir    Abundance results output directory (default: results_dir/abundance)

Example:
    bash batch_calculate_abundance_en.sh results_viral_short
    bash batch_calculate_abundance_en.sh results_viral_long abundance_output

Output:
    One abundance table file per sample: <sample>_abundance.tsv
    Summary file: all_samples_abundance_summary.tsv

EOF
    exit 1
}

# Check parameters
if [ $# -lt 1 ]; then
    usage
fi

RESULTS_DIR=$1
OUTPUT_DIR=${2:-"${RESULTS_DIR}/abundance"}

# Check if results directory exists
if [ ! -d "$RESULTS_DIR" ]; then
    echo "❌ Error: Results directory does not exist: $RESULTS_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "🧬 Batch Calculate Viral Abundance"
echo "=========================================="
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Find all Bracken output files
BRACKEN_DIR="${RESULTS_DIR}/bracken"
KRAKEN_DIR="${RESULTS_DIR}/kraken2"

if [ ! -d "$KRAKEN_DIR" ]; then
    echo "❌ Error: Kraken2 results directory not found: $KRAKEN_DIR"
    exit 1
fi

# Check if Bracken directory exists and is not empty
if [ ! -d "$BRACKEN_DIR" ] || [ -z "$(find "$BRACKEN_DIR" -name "*.tsv" 2>/dev/null)" ]; then
    echo "⚠️  Warning: Bracken results not found"
    echo ""
    echo "   Bracken may not have run or database is misconfigured."
    echo "   Will calculate abundance directly from Kraken2 results (similar to long-read method)."
    echo ""
    echo "   💡 This will still give you valid abundance values, just without Bracken statistical correction."
    echo ""
    
    # Fallback to long-read script
    if [ -f "batch_calculate_abundance_longread_en.sh" ]; then
        echo "   Using alternative method: batch_calculate_abundance_longread_en.sh"
        bash batch_calculate_abundance_longread_en.sh "$RESULTS_DIR" "$OUTPUT_DIR"
        exit $?
    else
        echo "❌ Error: batch_calculate_abundance_longread_en.sh not found either"
        echo "   Please ensure abundance calculation scripts are uploaded to the server"
        exit 1
    fi
fi

echo "✅ Found Bracken results, will use standard method for abundance calculation"

# Counters
processed=0
failed=0

# Find all possible Bracken output file formats (including subdirectories)
shopt -s nullglob  # If no matches, don't return pattern itself
BRACKEN_FILES=()

# Search in main directory
for file in "${BRACKEN_DIR}"/*_bracken*.tsv "${BRACKEN_DIR}"/*.bracken.tsv "${BRACKEN_DIR}"/*.bracken_species.tsv; do
    [ -f "$file" ] && BRACKEN_FILES+=("$file")
done

# Search in subdirectories (important! nf-core/taxprofiler may output to subdirectories)
for file in "${BRACKEN_DIR}"/*/*_bracken*.tsv "${BRACKEN_DIR}"/*/*.bracken.tsv "${BRACKEN_DIR}"/*/*.bracken_species.tsv; do
    [ -f "$file" ] && BRACKEN_FILES+=("$file")
done

# Also try finding any .tsv files in Bracken directory (in case naming is different)
if [ ${#BRACKEN_FILES[@]} -eq 0 ]; then
    while IFS= read -r file; do
        [ -f "$file" ] && BRACKEN_FILES+=("$file")
    done < <(find "$BRACKEN_DIR" -type f -name "*.tsv" 2>/dev/null | head -20)
fi

if [ ${#BRACKEN_FILES[@]} -eq 0 ]; then
    echo "⚠️  Warning: No Bracken output files found with expected patterns"
    echo "   Tried patterns: *_bracken*.tsv, *.bracken.tsv, *.bracken_species.tsv"
    echo "   Tried locations: ${BRACKEN_DIR}/ and subdirectories"
    echo ""
    echo "   Bracken directory contents:"
    ls -lh "$BRACKEN_DIR/" 2>/dev/null || echo "   Directory empty or doesn't exist"
    echo ""
    echo "   All .tsv files in Bracken directory:"
    find "$BRACKEN_DIR" -type f -name "*.tsv" 2>/dev/null | head -10 || echo "   No .tsv files found"
    echo ""
    echo "   Falling back to Kraken2-only method..."
    if [ -f "batch_calculate_abundance_longread_en.sh" ]; then
        bash batch_calculate_abundance_longread_en.sh "$RESULTS_DIR" "$OUTPUT_DIR"
        exit $?
    else
        echo "❌ Error: batch_calculate_abundance_longread_en.sh not found"
        exit 1
    fi
fi

echo "✅ Found ${#BRACKEN_FILES[@]} Bracken output file(s)"
echo ""

# Iterate through all found Bracken output files
for bracken_file in "${BRACKEN_FILES[@]}"; do
    if [ ! -f "$bracken_file" ]; then
        continue
    fi
    
    # Extract sample name (handle various possible file name formats)
    sample=$(basename "$bracken_file" | \
             sed 's/_bracken.*\.tsv$//' | \
             sed 's/\.bracken.*\.tsv$//' | \
             sed 's/_null_.*$//' | \
             sed 's/_.*_Viral_ref$//' | \
             sed 's/_Viral_ref$//')
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "📊 Processing sample: $sample"
    
    # Find corresponding Kraken2 report file (search in main directory and subdirectories)
    kraken_file=""
    
    # Try various naming formats in main directory
    for pattern in "${sample}.report" "${sample}_kraken2.report" "${sample}.kraken2.report.txt" "${sample}_kraken2_report.txt"; do
        if [ -f "${KRAKEN_DIR}/${pattern}" ]; then
            kraken_file="${KRAKEN_DIR}/${pattern}"
            break
        fi
    done
    
    # If not found, search in subdirectories
    if [ -z "$kraken_file" ]; then
        kraken_file=$(find "$KRAKEN_DIR" -type f \( -name "${sample}*.report" -o -name "${sample}*.kreport" -o -name "${sample}*.kraken2.report.txt" -o -name "${sample}*_kraken2_report.txt" \) | head -n 1)
    fi
    
    if [ ! -f "$kraken_file" ]; then
        echo "❌ Error: Kraken2 report file not found: $kraken_file"
        ((failed++))
        continue
    fi
    
    # Output file
    output_file="${OUTPUT_DIR}/${sample}_abundance.tsv"
    
    # Run calculation script
    echo "   Bracken: $bracken_file"
    echo "   Kraken:  $kraken_file"
    echo "   Output:  $output_file"
    echo "   Python:  $PYTHON_CMD"
    
    # Verify that the Python command is available
    if [ ! -f "$PYTHON_CMD" ] && ! command -v "$PYTHON_CMD" &> /dev/null; then
        echo "❌ Error: Python command not found: $PYTHON_CMD"
        ((failed++))
        continue
    fi
    
    "$PYTHON_CMD" calculate_abundance_en.py \
        --bracken "$bracken_file" \
        --kraken "$kraken_file" \
        --output "$output_file"
    
    if [ $? -eq 0 ]; then
        echo "✅ Complete: $sample"
        ((processed++))
    else
        echo "❌ Failed: $sample"
        ((failed++))
    fi
    
    echo ""
done

echo "=========================================="
echo "📈 Processing Complete"
echo "=========================================="
echo "Successfully processed: $processed samples"
echo "Failed: $failed samples"
echo "Output directory: $OUTPUT_DIR"

# If samples were successfully processed, generate summary table
if [ $processed -gt 0 ]; then
    echo ""
    echo "📋 Generating summary table..."
    
    summary_file="${OUTPUT_DIR}/all_samples_abundance_summary.tsv"
    
    # Create header
    echo -e "Sample\tSpecies\tTaxonomy_ID\tAssigned_Reads\tFraction\tRPM" > "$summary_file"
    
    # Merge all sample results
    for abundance_file in ${OUTPUT_DIR}/*_abundance.tsv; do
        if [ -f "$abundance_file" ]; then
            sample=$(basename "$abundance_file" | sed 's/_abundance\.tsv//')
            # Skip header, add sample column
            tail -n +2 "$abundance_file" | awk -v s="$sample" '{print s"\t"$0}' >> "$summary_file"
        fi
    done
    
    echo "✅ Summary table generated: $summary_file"
    
    # Generate TOP virus summary
    top_summary="${OUTPUT_DIR}/top_viruses_summary.tsv"
    echo ""
    echo "📊 Generating TOP virus summary (RPM >= 10)..."
    
    # Use awk to extract viruses with RPM >= 10 and sort
    awk -F'\t' 'NR==1 || $6 >= 10' "$summary_file" | \
        sort -t$'\t' -k6 -nr > "$top_summary"
    
    echo "✅ TOP virus summary: $top_summary"
fi

echo ""
echo "🎉 All tasks complete!"

