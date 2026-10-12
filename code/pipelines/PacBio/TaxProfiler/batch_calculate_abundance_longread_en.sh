#!/bin/bash
# Batch calculate long-read data abundance (based on Kraken2 report)

# === Clean a polluted Python environment ===
# The login environment exports PYTHONPATH/PYTHONHOME pointing at the old
# user-py39 site-packages (e.g. PYTHONPATH=:/home/$USER/user-py39/lib/
# python3.9/site-packages). These leak into the conda Python selected below and
# break interpreter startup ("init_fs_encoding: failed to get the Python codec
# of the filesystem encoding"), so `import pandas, numpy` fails even after they
# are installed. Unset them before touching any python/pip/conda.
unset PYTHONPATH PYTHONHOME

# === Detect and use correct Python interpreter ===
# Prefer the Python from the conda environment so that pandas and numpy are available

# Build an ORDERED list of candidate interpreters. We deliberately do NOT commit
# to one by mere file existence: a conda env can contain a python that cannot
# even start (e.g. an env built for Nextflow/Java whose stdlib is incomplete —
# the "init_fs_encoding / can't find encoding" crash). The real choice is made
# further down by testing which candidate can actually `import pandas, numpy`.
CANDIDATES=()
[ -n "${PYTHON_CMD:-}" ]    && CANDIDATES+=("$PYTHON_CMD")            # parent-provided
[ -n "${CONDA_PREFIX:-}" ]  && CANDIDATES+=("$CONDA_PREFIX/bin/python")
CANDIDATES+=(
    "$HOME/.conda/envs/nextflow_env/bin/python"
    "/scratch/$USER/conda/envs/nextflow_env/bin/python"
    "/apps/eb/Miniforge3/24.11.3-0/envs/nextflow_env/bin/python"
    "/opt/conda/envs/nextflow_env/bin/python"
)
command -v python3 &> /dev/null && CANDIDATES+=("$(command -v python3)")
command -v python  &> /dev/null && CANDIDATES+=("$(command -v python)")

# === Select the FIRST candidate that actually works AND has pandas + numpy ===
# The import test doubles as an interpreter health check: a broken python (one
# that cannot load its own stdlib) fails `import pandas, numpy` exactly like a
# python that is missing the package. A single test therefore filters out BOTH
# failure modes — which is why "just installing pandas" did not help before: the
# selected interpreter could not start at all.
#
# NOTE: automatic pip/conda installation was removed on purpose. Installing
# packages inside a batch job — against an interpreter that may be broken — only
# masks the real problem and wastes wall-clock. If nothing works we fail loudly
# with the list of interpreters tried.
PYTHON_CMD=""
echo "🔍 Selecting a working Python with pandas + numpy ..."
for cand in "${CANDIDATES[@]}"; do
    [ -n "$cand" ] || continue
    { [ -f "$cand" ] || command -v "$cand" &> /dev/null; } || continue
    if "$cand" -c "import pandas, numpy" 2>/dev/null; then
        PYTHON_CMD="$cand"
        echo "✅ Using: $PYTHON_CMD"
        "$PYTHON_CMD" -c "import sys, pandas, numpy; print('   executable:', sys.executable); print(f'   pandas {pandas.__version__}, numpy {numpy.__version__}')"
        break
    else
        echo "   ✗ skip (python broken or no pandas/numpy): $cand"
    fi
done

if [ -z "$PYTHON_CMD" ]; then
    echo ""
    echo "❌ No Python interpreter with a WORKING pandas + numpy was found."
    echo "   Candidates tried (in order):"
    for cand in "${CANDIDATES[@]}"; do echo "     - $cand"; done
    echo ""
    echo "   Fix: install pandas/numpy into a working conda env, then confirm it"
    echo "   actually runs (not just that files exist):"
    echo "     conda activate nextflow_env"
    echo "     conda install -c conda-forge pandas numpy -y"
    echo "     \"\$(command -v python)\" -c 'import pandas, numpy; print(pandas.__version__)'"
    echo ""
    echo "   If one of the envs above prints the fatal 'init_fs_encoding /"
    echo "   can''t find encoding' error, that env's python is broken/incomplete"
    echo "   — recreate it or use a different one; do NOT point this script at it."
    exit 1
fi
echo ""

# Usage instructions
usage() {
    cat << EOF
Usage: bash batch_calculate_abundance_longread_en.sh <results_dir> [output_dir]

Parameters:
    results_dir   TaxProfiler output results directory (e.g. results_viral_long)
    output_dir    Abundance results output directory (default: results_dir/abundance)

Description:
    Long-read data doesn't use Bracken, extracts abundance directly from Kraken2 report.

Example:
    bash batch_calculate_abundance_longread_en.sh results_viral_long

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
echo "🧬 Batch Calculate Long-read Viral Abundance (Kraken2)"
echo "=========================================="
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Find all Kraken2 report files
KRAKEN_DIR="${RESULTS_DIR}/kraken2"

if [ ! -d "$KRAKEN_DIR" ]; then
    echo "❌ Error: Kraken2 results directory not found: $KRAKEN_DIR"
    exit 1
fi

# Counters
processed=0
failed=0

# Find all possible Kraken2 report file formats (including subdirectories)
shopt -s nullglob  # If no matches, don't return pattern itself
REPORT_FILES=()

# Search in main directory
REPORT_FILES+=("${KRAKEN_DIR}"/*.report)
REPORT_FILES+=("${KRAKEN_DIR}"/*.kreport)
REPORT_FILES+=("${KRAKEN_DIR}"/*.kraken2.report.txt)
REPORT_FILES+=("${KRAKEN_DIR}"/*_kraken2_report.txt)

# Search in subdirectories (important! nf-core/taxprofiler may output to subdirectories)
REPORT_FILES+=("${KRAKEN_DIR}"/*/*.report)
REPORT_FILES+=("${KRAKEN_DIR}"/*/*.kreport)
REPORT_FILES+=("${KRAKEN_DIR}"/*/*.kraken2.report.txt)
REPORT_FILES+=("${KRAKEN_DIR}"/*/*_kraken2_report.txt)

if [ ${#REPORT_FILES[@]} -eq 0 ]; then
    echo "⚠️  Kraken2 report files not found"
    echo "   Tried file formats: *.report, *.kreport, *.kraken2.report.txt"
    echo "   Tried locations: ${KRAKEN_DIR}/ and subdirectories"
    echo ""
    echo "   Kraken2 directory contents:"
    ls -lh "$KRAKEN_DIR/" 2>/dev/null || echo "   Directory empty or doesn't exist"
    echo ""
    echo "   Subdirectory contents:"
    find "$KRAKEN_DIR" -type f -name "*.txt" 2>/dev/null | head -5
    exit 1
fi

echo "✅ Found ${#REPORT_FILES[@]} Kraken2 report file(s)"
echo ""

# Iterate through all found Kraken2 report files
for kraken_file in "${REPORT_FILES[@]}"; do
    # Extract sample name (handle various possible file name formats)
    sample=$(basename "$kraken_file" | \
             sed 's/\.kraken2\.report\.txt$//' | \
             sed 's/_kraken2_report\.txt$//' | \
             sed 's/\.report$//' | \
             sed 's/\.kreport$//' | \
             sed 's/_null_.*$//' | \
             sed 's/_.*_Viral_ref$//' | \
             sed 's/_kraken2$//')
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "📊 Processing sample: $sample"
    echo "   Kraken2: $kraken_file"
    
    # Output file
    output_file="${OUTPUT_DIR}/${sample}_abundance.tsv"
    echo "   Output:  $output_file"
    
    # Run calculation script
    "$PYTHON_CMD" calculate_abundance_longread_en.py \
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

