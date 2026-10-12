#!/bin/bash
#SBATCH --job-name=CLARK_sam-2_short
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=CLARK_708_short_sam-2_%j.out
#SBATCH --error=CLARK_708_short_sam-2_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# Edited: 21, 22, 23, 24
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/CLARK_DB"}
# CLARK_DIR: CLARK software installation directory, should contain classify_metagenome.sh, set_targets.sh, estimate_abundance.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/$USER/Meta-genome-data-analysis/Apptainer/CLARK"}
INPUT_R1=${INPUT_R1:-"/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/illumina-708-genomes/sample_708_2_01.fq.gz"}  # Paired-end R1 file (gzipped)
INPUT_R2=${INPUT_R2:-"/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/illumina-708-genomes/sample_708_2_02.fq.gz"}  # Paired-end R2 file (gzipped)
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"virus_results_708_2"}  # Output file prefix
OUTPUT_DIR=${OUTPUT_DIR:-"results_short-708_2"}  # Output directory
K_MER=${K_MER:-20}  # k-mer size (recommended 20-21 for virus detection, higher sensitivity)
MODE=${MODE:-0}  # Execution mode: 0=full, 1=default, 2=express. MUST be 0 to get the Gamma column.
GAMMA_MIN=${GAMMA_MIN:-0.2}  # post-processing: keep only reads with gamma >= this value before computing abundance (removes bacterial noise, improves precision).
                             # Requires MODE=0 (full) for the gamma column. Short reads may need the best threshold to be rescanned.
THREADS=$SLURM_CPUS_PER_TASK  # Use CPU count allocated by SLURM

# Abundance estimation parameters (optional)
ENABLE_ABUNDANCE=${ENABLE_ABUNDANCE:-"true"}  # whether to run abundance estimation after classification
CONFIDENCE_THRESHOLD=${CONFIDENCE_THRESHOLD:-""}  # confidence threshold, for example: 0.80
GAMMA_THRESHOLD=${GAMMA_THRESHOLD:-""}  # Gamma threshold, for example: 0.03
ABUNDANCE_THRESHOLD=${ABUNDANCE_THRESHOLD:-""}  # abundance threshold (%), for example: 2
HIGH_CONFIDENCE=${HIGH_CONFIDENCE:-"false"}  # whether to use only high-confidence assignments
MPA_FORMAT=${MPA_FORMAT:-"false"}  # whether to output MetaPhlAn format

# ============================================
# Start Classification
# ============================================
echo "=========================================="
echo "CLARK Virus Classification Job (Paired-end)"
echo "Time: $(date)"
echo "Input file R1: $INPUT_R1"
echo "Input file R2: $INPUT_R2"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Database directory: $CLARK_DB_DIR"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "Threads: $THREADS"
echo "Abundance estimation: $ENABLE_ABUNDANCE"
if [ "$ENABLE_ABUNDANCE" = "true" ]; then
    [ -n "$CONFIDENCE_THRESHOLD" ] && echo "  Confidence threshold: $CONFIDENCE_THRESHOLD"
    [ -n "$GAMMA_THRESHOLD" ] && echo "  Gamma threshold: $GAMMA_THRESHOLD"
    [ -n "$ABUNDANCE_THRESHOLD" ] && echo "  Abundance threshold: $ABUNDANCE_THRESHOLD%"
    [ "$HIGH_CONFIDENCE" = "true" ] && echo "  High confidence mode: enabled"
    [ "$MPA_FORMAT" = "true" ] && echo "  MetaPhlAn format: enabled"
fi
echo "=========================================="

# Check input files
if [ ! -f "$INPUT_R1" ]; then
    echo "Error: R1 file does not exist: $INPUT_R1"
    exit 1
fi

if [ ! -f "$INPUT_R2" ]; then
    echo "Error: R2 file does not exist: $INPUT_R2"
    exit 1
fi

# Check and validate the database
# CLARK database structure: CLARK_db/viruses_*/db_central_k*_*.tsk.*
# Or:CLARK_db/Viruses/viruses_*/db_central_k*_*.tsk.*
validate_database() {
    local db_dir=$1
    local found_db_path=""
    local found_db_files=""
    
    echo "Validating CLARK database..."
    
    # Check whether the database directory exists
    if [ ! -d "$db_dir" ]; then
        echo "ERROR: Database directory does not exist: $db_dir"
        return 1
    fi
    
    echo "  Database root directory: $db_dir"
    
    # Use the find command to search recursively for the database files (the most reliable method)
    # CLARK database file naming format: db_central_k{K_MER}_*.tsk.*
    echo "  Searching for database files with k-mer size $K_MER..."
    found_db_files=$(find "$db_dir" -type f -name "db_central_k${K_MER}_*.tsk.*" 2>/dev/null | head -1)
    
    if [ -n "$found_db_files" ]; then
        found_db_path=$(dirname "$found_db_files")
    fi
    
    # Check whether database files were found
    if [ -z "$found_db_files" ]; then
        echo "ERROR: Cannot find database files for k-mer size $K_MER"
        echo ""
        echo "Searched in: $db_dir"
        echo "Looking for pattern: db_central_k${K_MER}_*.tsk.*"
        echo ""
        
        # Show the database directory structure
        echo "Database directory structure:"
        ls -la "$db_dir" 2>/dev/null | head -20 || echo "  (Cannot list directory contents)"
        echo ""
        
        # Show the available database files
        echo "Available database files in $db_dir:"
        local available_files=$(find "$db_dir" -type f -name "db_central_k*.tsk.*" 2>/dev/null | head -10)
        if [ -n "$available_files" ]; then
            echo "$available_files" | while read -r f; do
                echo "  - $f"
            done
        else
            echo "  (No database files found)"
        fi
        echo ""
        
        # Show the available k-mer sizes
        echo "Available k-mer sizes:"
        local kmer_sizes=$(find "$db_dir" -type f -name "db_central_k*.tsk.*" 2>/dev/null | sed 's/.*db_central_k\([0-9]*\)_.*/\1/' | sort -u)
        if [ -n "$kmer_sizes" ]; then
            echo "$kmer_sizes" | while read -r k; do
                echo "  - k=$k"
            done
        else
            echo "  (No k-mer sizes found)"
        fi
        echo ""
        
        echo "Please ensure:"
        echo "  1. Database is built for k-mer size $K_MER"
        echo "  2. Database files are in the correct location"
        echo "  3. Run slurm_build_virus_db.sh if database is not built"
        return 1
    fi
    
    # Verify the database files found
    echo "  Database files found!"
    echo "  Database path: $found_db_path"
    echo "  Sample database file: $(basename "$found_db_files")"
    
    # Check the database file size
    local db_size=$(du -sh "$found_db_path" 2>/dev/null | awk '{print $1}')
    if [ -n "$db_size" ]; then
        echo "  Database directory size: $db_size"
    fi
    
    # List all matching database files (for verification)
    local file_count=$(find "$found_db_path" -type f -name "db_central_k${K_MER}_*.tsk.*" 2>/dev/null | wc -l)
    if [ "$file_count" -gt 0 ]; then
        echo "  Found $file_count database file(s) for k-mer size $K_MER"
    fi
    
    return 0
}

if ! validate_database "$CLARK_DB_DIR"; then
    # db_central index not built yet. If CLARK targets are defined (set_targets.sh
    # was run -> .viruses/.bacteria + targets.txt exist), CLARK will BUILD the
    # database automatically on the first classify run. Proceed instead of aborting.
    if [ -f "$CLARK_DB_DIR/targets.txt" ] || ls "$CLARK_DB_DIR"/.viruses >/dev/null 2>&1; then
        echo "=========================================="
        echo "NOTE: db_central (k=$K_MER) not found, but CLARK targets ARE defined."
        echo "      CLARK will BUILD the database on this first classify run."
        echo "      WARNING: building a multi-domain DB (bacteria+viruses+...) at"
        echo "      k=$K_MER is RAM/disk/time heavy. Ensure --mem and --time are large"
        echo "      enough, or the job may OOM / time out mid-build."
        echo "=========================================="
    else
        echo "ERROR: no built DB and no targets defined in $CLARK_DB_DIR."
        echo "       Run set_targets.sh first, e.g.:"
        echo "       set_targets.sh \$DIR_DB bacteria viruses --species"
        exit 1
    fi
fi

# Check the CLARK installation
if [ ! -f "$CLARK_DIR/classify_metagenome.sh" ]; then
    echo "Error: Cannot find classify_metagenome.sh script!"
    echo "Please set CLARK_DIR environment variable"
    exit 1
fi

# Check the CLARK executable
if [ ! -f "$CLARK_DIR/exe/CLARK" ]; then
    echo "WARNING: Cannot find CLARK executable at $CLARK_DIR/exe/CLARK"
    echo "CLARK may not be properly installed"
else
    # Check the executable permissions
    if [ ! -x "$CLARK_DIR/exe/CLARK" ]; then
        echo "WARNING: CLARK executable is not executable: $CLARK_DIR/exe/CLARK"
        echo "Attempting to fix permissions..."
        chmod +x "$CLARK_DIR/exe/CLARK" 2>/dev/null || {
            echo "ERROR: Cannot make CLARK executable. Please check permissions."
            exit 1
        }
    fi
    echo "CLARK executable found and verified: $CLARK_DIR/exe/CLARK"
fi

# If abundance estimation is enabled, check the related tools
if [ "$ENABLE_ABUNDANCE" = "true" ]; then
    if [ -f "$CLARK_DIR/estimate_abundance.sh" ]; then
        echo "estimate_abundance.sh found - will use it for abundance estimation"
    elif command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
        echo "estimate_abundance.sh not found, but Python is available - will use manual calculation"
    else
        echo "WARNING: estimate_abundance.sh not found and Python is not available"
        echo "Abundance estimation may fail. Please install Python 3 or ensure estimate_abundance.sh exists."
    fi
fi

# Convert OUTPUT_DIR to absolute path if it's relative
if [[ "$OUTPUT_DIR" != /* ]]; then
    # Relative path - make it relative to submission directory
    OUTPUT_DIR_ABS="$SLURM_SUBMIT_DIR/$OUTPUT_DIR"
else
    OUTPUT_DIR_ABS="$OUTPUT_DIR"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR_ABS"

# Set TMPDIR to a location with sufficient space (use scratch directory instead of /tmp)
# CLARK needs to decompress gzipped files, which requires temporary space
# Force TMPDIR to use scratch directory to avoid /tmp space issues
# Note: This overrides any system default TMPDIR setting
TMPDIR="/scratch/$USER/tmp"
export TMPDIR
mkdir -p "$TMPDIR"
echo "Temporary directory (TMPDIR) forced to: $TMPDIR"
# Check available space in TMPDIR (in GB)
AVAILABLE_SPACE=$(df -BG "$TMPDIR" 2>/dev/null | tail -1 | awk '{print $4}' | sed 's/G//')
if [ -n "$AVAILABLE_SPACE" ] && [ "$AVAILABLE_SPACE" -lt 10 ]; then
    echo "WARNING: Low disk space in $TMPDIR: ${AVAILABLE_SPACE}GB available"
    echo "CLARK may fail when decompressing large gzipped files"
fi

# Make sure the CLARK database environment variable is set
# CLARK 's classify_metagenome.sh uses the CLARK_DB environment variable (not CLARK_DB_DIR)
# but some versions may use CLARK_DB_DIR, so set both
export CLARK_DB="$CLARK_DB_DIR"
export CLARK_DB_DIR="$CLARK_DB_DIR"

echo "CLARK database environment variables:"
echo "  CLARK_DB=$CLARK_DB"
echo "  CLARK_DB_DIR=$CLARK_DB_DIR"

# Enter CLARK directory
cd "$CLARK_DIR" || exit 1

# Build classification command (paired-end data)
# Note: classify_metagenome.sh uses -P for paired-end mode
echo "Start time: $(date)"
echo "----------------------------------------"

# Set memory limits and system resource limits to prevent out-of-memory segmentation faults
# Set the virtual memory limit (about 90% of the allocated memory, leaving a buffer)
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    # SLURM_MEM_PER_NODE unit is MB, convert to KB
    MEM_LIMIT_KB=$((SLURM_MEM_PER_NODE * 1024 * 9 / 10))
    ulimit -v $MEM_LIMIT_KB 2>/dev/null || echo "WARNING: Could not set virtual memory limit"
    echo "Set virtual memory limit to ${MEM_LIMIT_KB}KB (90% of allocated ${SLURM_MEM_PER_NODE}MB)"
fi

# Set the maximum file size limit (prevents creating oversized temporary files)
ulimit -f unlimited 2>/dev/null || echo "WARNING: Could not set file size limit"

# Check the available memory
if command -v free >/dev/null 2>&1; then
    AVAILABLE_MEM=$(free -g | awk '/^Mem:/ {print $7}')
    echo "Available system memory: ${AVAILABLE_MEM}GB"
    if [ -n "$AVAILABLE_MEM" ] && [ "$AVAILABLE_MEM" -lt 20 ]; then
        echo "WARNING: Low available memory (${AVAILABLE_MEM}GB). Consider reducing thread count."
    fi
fi

# Build the command array for safe execution
# Note: the -O option accepts multiple files for paired-end data
# Compute a safe number of threads (use 50% of the available threads to avoid memory problems)
# CLARK has very high memory requirements, so the number of threads is limited to prevent segmentation faults
# From experience, 4-6 threads are usually more stable for large files and large databases
SAFE_THREADS=$((THREADS / 2))
if [ $SAFE_THREADS -lt 1 ]; then
    SAFE_THREADS=1
fi
# Limit further to 4 threads for stability (adjust according to memory)
# If it failed before, you can manually set fewer threads
if [ $SAFE_THREADS -gt 4 ]; then
    SAFE_THREADS=4
fi

# If the environment variable FORCE_THREADS is set, use that value
if [ -n "$FORCE_THREADS" ] && [ "$FORCE_THREADS" -gt 0 ]; then
    SAFE_THREADS=$FORCE_THREADS
    echo "Using forced thread count from FORCE_THREADS environment variable"
fi

echo "Using $SAFE_THREADS threads (reduced from $THREADS to prevent memory issues)"

# Handle compressed files: CLARK's classify_metagenome.sh may hardcode /tmp for decompression
# To avoid /tmp space issues, we'll decompress files to output directory first
UNCOMPRESSED_R1="$INPUT_R1"
UNCOMPRESSED_R2="$INPUT_R2"
CLEANUP_FILES=()

# Function to clean up decompressed files on exit
cleanup_decompressed() {
    if [ ${#CLEANUP_FILES[@]} -gt 0 ]; then
        echo "Cleaning up decompressed temporary files..."
        for file in "${CLEANUP_FILES[@]}"; do
            if [ -f "$file" ]; then
                rm -f "$file"
                echo "  Removed: $file"
            fi
        done
    fi
}

# Set trap to cleanup on exit (normal or error)
trap cleanup_decompressed EXIT

# Verify the integrity of the input file
validate_input_file() {
    local file=$1
    local file_type=$2
    
    echo "Validating $file_type file: $file"
    
    # Check whether the file is readable
    if [ ! -r "$file" ]; then
        echo "ERROR: Cannot read $file_type file: $file"
        return 1
    fi
    
    # Check the file size
    local file_size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo "0")
    if [ "$file_size" -eq 0 ]; then
        echo "ERROR: $file_type file is empty: $file"
        return 1
    fi
    # Format the file size display (compatible with different systems)
    local size_display
    if command -v numfmt >/dev/null 2>&1; then
        size_display=$(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null)
    else
        # Format manually (convert to GB/MB/KB)
        if [ $file_size -gt 1073741824 ]; then
            size_display=$(awk "BEGIN {printf \"%.2fGB\", $file_size/1073741824}")
        elif [ $file_size -gt 1048576 ]; then
            size_display=$(awk "BEGIN {printf \"%.2fMB\", $file_size/1048576}")
        elif [ $file_size -gt 1024 ]; then
            size_display=$(awk "BEGIN {printf \"%.2fKB\", $file_size/1024}")
        else
            size_display="${file_size} bytes"
        fi
    fi
    echo "  File size: $size_display"
    
    # If it is a compressed file, verify the integrity of the archive
    if [[ "$file" == *.gz ]]; then
        echo "  Checking gzip integrity..."
        if ! gzip -t "$file" 2>/dev/null; then
            echo "ERROR: $file_type file is corrupted (gzip test failed): $file"
            return 1
        fi
        echo "  Gzip integrity check passed"
    # If it is a FASTQ file, do a basic format check
    elif [[ "$file" == *.fastq ]] || [[ "$file" == *.fq ]]; then
        echo "  Checking FASTQ format..."
        local line_count=$(head -n 4 "$file" 2>/dev/null | wc -l)
        if [ "$line_count" -lt 4 ]; then
            echo "WARNING: $file_type file may not be in valid FASTQ format (less than 4 lines in header)"
        else
            local first_line=$(head -n 1 "$file" 2>/dev/null)
            if [[ "$first_line" != @* ]]; then
                echo "WARNING: $file_type file may not be in valid FASTQ format (first line doesn't start with @)"
            else
                echo "  FASTQ format check passed"
            fi
        fi
    fi
    
    return 0
}

# Verify the input file
if ! validate_input_file "$INPUT_R1" "R1"; then
    echo "ERROR: R1 file validation failed"
    exit 1
fi

if ! validate_input_file "$INPUT_R2" "R2"; then
    echo "ERROR: R2 file validation failed"
    exit 1
fi

if [[ "$INPUT_R1" == *.gz ]] || [[ "$INPUT_R2" == *.gz ]]; then
    echo "Detected gzip compressed files. Decompressing to avoid /tmp space issues..."
    
    # Decompress R1 (if needed)
    if [[ "$INPUT_R1" == *.gz ]]; then
        UNCOMPRESSED_R1="$OUTPUT_DIR_ABS/$(basename "$INPUT_R1" .gz)"
        echo "Decompressing $INPUT_R1 to $UNCOMPRESSED_R1..."
        if ! gunzip -c "$INPUT_R1" > "$UNCOMPRESSED_R1"; then
            echo "Error: Failed to decompress $INPUT_R1"
            exit 1
        fi
        # Verify the decompressed files
        if ! validate_input_file "$UNCOMPRESSED_R1" "decompressed R1"; then
            echo "Error: Decompressed R1 file validation failed"
            rm -f "$UNCOMPRESSED_R1"
            exit 1
        fi
        CLEANUP_FILES+=("$UNCOMPRESSED_R1")
    fi
    
    # Decompress R2 (if needed)
    if [[ "$INPUT_R2" == *.gz ]]; then
        UNCOMPRESSED_R2="$OUTPUT_DIR_ABS/$(basename "$INPUT_R2" .gz)"
        echo "Decompressing $INPUT_R2 to $UNCOMPRESSED_R2..."
        if ! gunzip -c "$INPUT_R2" > "$UNCOMPRESSED_R2"; then
            echo "Error: Failed to decompress $INPUT_R2"
            # Clean up the decompressed R1 (if it was decompressed)
            [ -f "$UNCOMPRESSED_R1" ] && rm -f "$UNCOMPRESSED_R1"
            exit 1
        fi
        # Verify the decompressed files
        if ! validate_input_file "$UNCOMPRESSED_R2" "decompressed R2"; then
            echo "Error: Decompressed R2 file validation failed"
            rm -f "$UNCOMPRESSED_R2"
            [ -f "$UNCOMPRESSED_R1" ] && rm -f "$UNCOMPRESSED_R1"
            exit 1
        fi
        CLEANUP_FILES+=("$UNCOMPRESSED_R2")
    fi
    
    echo "Decompression complete. Using uncompressed files for CLARK."
else
    echo "Input files are not compressed, using directly."
fi

# Function that runs the classification, with a retry mechanism
run_classification() {
    local threads=$1
    local attempt_num=$2
    
    echo "----------------------------------------"
    echo "Classification attempt #$attempt_num with $threads threads"
    echo "----------------------------------------"
    
    # Build the classification command (after the final file path has been determined)
    local classify_cmd=(
        "./classify_metagenome.sh"
        "-P"           # paired-end mode
        "$UNCOMPRESSED_R1"
        "$UNCOMPRESSED_R2"
        "-R"
        "$OUTPUT_DIR_ABS/$OUTPUT_PREFIX"
        "-m"
        "$MODE"
        "-k"
        "$K_MER"
        "-n"
        "$threads"
    )
    
    echo "CLARK command will use:"
    echo "  R1: $UNCOMPRESSED_R1"
    echo "  R2: $UNCOMPRESSED_R2"
    echo "  Output prefix: $OUTPUT_DIR_ABS/$OUTPUT_PREFIX"
    echo "  Threads: $threads"
    echo "  Mode: $MODE"
    echo "  k-mer size: $K_MER"
    
    # Run the classification
    # Note: classify_metagenome.sh may print a "Failed to recognize option" warning for the R2 file
    # This is a false alarm - the script still processes the paired-end data correctly
    # Capture stdout and stderr and filter warnings, but keep the exit code
    # Use a temporary file to capture the output and check for errors
    local temp_output=$(mktemp)
    "${classify_cmd[@]}" 2>&1 | grep -v "^Failed to recognize option:" | tee "$temp_output"
    local exit_code=${PIPESTATUS[0]}
    
    # Check for the various error conditions in the output and record the error type
    local error_type=""
    if grep -qi "segmentation fault\|core dumped" "$temp_output"; then
        exit_code=139  # SIGSEGV 's standard exit code
        error_type="segfault"
        echo "WARNING: Segmentation fault detected in CLARK execution!"
    elif grep -qi "No space left on device\|Failed to copy and uncompress\|Failed to.*uncompress" "$temp_output"; then
        exit_code=1  # general error code
        error_type="disk_space"
        echo "ERROR: Disk space issue detected - CLARK failed to decompress input files!"
        echo "The temporary directory ($TMPDIR) may be full or have insufficient space."
    elif grep -qi "Failed to\|Error\|error" "$temp_output"; then
        # Check whether it is a real error (not only the warnings we filter out)
        if ! grep -qi "Failed to recognize option" "$temp_output"; then
            # If the exit code is still 0 but errors were seen, set it to 1
            if [ $exit_code -eq 0 ]; then
                exit_code=1
                error_type="general"
                echo "WARNING: Error messages detected in CLARK output, but exit code was 0"
            fi
        fi
    fi
    
    # Clean up temporary files
    rm -f "$temp_output"
    
    # Return the exit code and error type (through global variables)
    CLASSIFY_EXIT_CODE=$exit_code
    ERROR_TYPE=$error_type
    
    return $exit_code
}

# Run the classification, with automatic retry support
MAX_RETRIES=2
CURRENT_THREADS=$SAFE_THREADS
ATTEMPT=1
CLASSIFY_EXIT_CODE=0
ERROR_TYPE=""

while [ $ATTEMPT -le $((MAX_RETRIES + 1)) ]; do
    if run_classification $CURRENT_THREADS $ATTEMPT; then
        # Completed successfully
        break
    elif [ "$ERROR_TYPE" = "segfault" ] && [ $ATTEMPT -le $MAX_RETRIES ]; then
        # Segmentation fault with retries left: retry with fewer threads
        echo "----------------------------------------"
        echo "Segmentation fault occurred. Retrying with fewer threads..."
        CURRENT_THREADS=$((CURRENT_THREADS / 2))
        if [ $CURRENT_THREADS -lt 1 ]; then
            CURRENT_THREADS=1
        fi
        echo "Reducing threads to $CURRENT_THREADS for retry attempt #$((ATTEMPT + 1))"
        ATTEMPT=$((ATTEMPT + 1))
        
        # Never fall back to express mode automatically (MODE=2)!
        # express mode outputs 3 columns without a gamma column, so gamma filtering cannot be done (precision cannot be improved).
        # We need the gamma of full mode (MODE=0), so keep MODE unchanged and only retry with fewer threads.
        # If full mode keeps segfaulting because of insufficient memory, raise --mem or export FORCE_THREADS=2,
        # rather than a degraded mode.
        if [ "$MODE" != "0" ]; then
            echo "WARNING: MODE is not 0 (full). Gamma column will be ABSENT -> cannot gamma-filter."
        fi
    else
        # Other errors, or retries exhausted: exit the loop
        break
    fi
done

# Note: Decompressed files will be cleaned up by the EXIT trap

echo "----------------------------------------"
echo "End time: $(date)"

if [ $CLASSIFY_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Virus classification completed successfully!"
    echo "Result file: $OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.csv"
    echo "=========================================="
    
    # Show information about the classification result file
    if [ -f "$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.csv" ]; then
        echo "Classification result file statistics:"
        wc -l "$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.csv"
    fi

    # ============================================
    # Gamma -filtered high-precision viral abundance table (requires the 8-column output of MODE=0)
    # ============================================
    RESULT_CSV_GF="$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.csv"
    if [ -f "$RESULT_CSV_GF" ] && [ -f "$CLARK_DIR/estimate_abundance.sh" ]; then
        NCOL=$(head -1 "$RESULT_CSV_GF" | awk -F, '{print NF}')
        if [ "$NCOL" -lt 8 ]; then
            echo ""
            echo "WARNING: result has only $NCOL columns (no Gamma column) -> not full-mode output, skipping gamma filtering."
            echo "         Please confirm MODE=0 and rerun; otherwise high-precision filtering cannot be done."
        else
            echo ""
            echo "=========================================="
            echo "Gamma filtered abundance estimate (gamma >= ${GAMMA_MIN})"
            echo "=========================================="
            ABUND_G="$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.g${GAMMA_MIN}.assignments.csv"
            # 8columns (full): Object_ID,Length,Gamma,1st_assignment,... -> take columns 1,2,4 and gamma>=threshold
            awk -F, -v G="$GAMMA_MIN" \
                'NR==1{print "Object_ID,Length,Assignment";next} ($3+0)>=G && $4!="NA"{print $1","$2","$4}' \
                "$RESULT_CSV_GF" > "$ABUND_G"
            KEPT=$(( $(wc -l < "$ABUND_G") - 1 ))
            echo "  reads passing gamma>=${GAMMA_MIN}: $KEPT"
            # Note: -D must come before -F (-F takes a variable number of arguments and would swallow the following -D)
            ( cd "$CLARK_DIR" && ./estimate_abundance.sh -D "$CLARK_DB_DIR" -F "$ABUND_G" ) \
                > "$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.abundance.g${GAMMA_MIN}.csv" \
                2> "$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.abundance.g${GAMMA_MIN}.log"
            echo "Wrote: $OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.abundance.g${GAMMA_MIN}.csv  <-- use this one to compute viral F1"
        fi
    fi

    # Run the abundance estimation (if enabled)
    if [ "$ENABLE_ABUNDANCE" = "true" ]; then
        echo ""
        echo "=========================================="
        echo "Starting abundance estimation..."
        echo "=========================================="
        
        # Initialize the variables related to abundance estimation
        ABUNDANCE_EXIT_CODE=0
        ABUNDANCE_OUTPUT_FILE=""
        
        # Check whether the classification result file exists
        CLASSIFICATION_CSV="$OUTPUT_DIR_ABS/${OUTPUT_PREFIX}.csv"
        if [ ! -f "$CLASSIFICATION_CSV" ]; then
            echo "WARNING: Classification result file not found: $CLASSIFICATION_CSV"
            echo "Skipping abundance estimation."
        else
            # Function that runs the abundance estimation
            # Two approaches are supported: 1) use estimate_abundance.sh (if present); 2) compute manually from the CSV
            run_abundance_estimation() {
                local input_csv=$1
                local input_base=$(basename "$input_csv" .csv)
                local abundance_output_file="$OUTPUT_DIR_ABS/${input_base}_abundance.csv"
                
                echo "Abundance estimation parameters:"
                echo "  Input file: $input_csv"
                [ -n "$CONFIDENCE_THRESHOLD" ] && echo "  Confidence threshold: $CONFIDENCE_THRESHOLD"
                [ -n "$GAMMA_THRESHOLD" ] && echo "  Gamma threshold: $GAMMA_THRESHOLD"
                [ -n "$ABUNDANCE_THRESHOLD" ] && echo "  Abundance threshold: $ABUNDANCE_THRESHOLD%"
                [ "$HIGH_CONFIDENCE" = "true" ] && echo "  High confidence mode: enabled"
                echo ""
                
                # Method 1: try estimate_abundance.sh (if present)
                if [ -f "$CLARK_DIR/estimate_abundance.sh" ]; then
                    echo "Using estimate_abundance.sh for abundance estimation..."
                    echo "Start time: $(date)"
                    echo "----------------------------------------"
                    
                    # Build the abundance estimation command
                    local abundance_cmd=(
                        "./estimate_abundance.sh"
                        "-F"
                        "$input_csv"
                        "-D"
                        "$CLARK_DB_DIR"
                    )
                    
                    # Add the filtering options
                    if [ "$HIGH_CONFIDENCE" = "true" ]; then
                        abundance_cmd+=("--highconfidence")
                    fi
                    
                    if [ -n "$CONFIDENCE_THRESHOLD" ]; then
                        abundance_cmd+=("-c" "$CONFIDENCE_THRESHOLD")
                    fi
                    
                    if [ -n "$GAMMA_THRESHOLD" ]; then
                        abundance_cmd+=("-g" "$GAMMA_THRESHOLD")
                    fi
                    
                    if [ -n "$ABUNDANCE_THRESHOLD" ]; then
                        abundance_cmd+=("-a" "$ABUNDANCE_THRESHOLD")
                    fi
                    
                    if [ "$MPA_FORMAT" = "true" ]; then
                        abundance_cmd+=("--mpa")
                    fi
                    
                    echo "Executing: ${abundance_cmd[*]}"
                    
                    # Run the abundance estimation
                    local temp_abundance_output=$(mktemp)
                    "${abundance_cmd[@]}" 2>&1 | tee "$temp_abundance_output"
                    local abundance_exit_code=${PIPESTATUS[0]}
                    
                    echo "----------------------------------------"
                    echo "End time: $(date)"
                    
                    # Check for errors
                    if grep -qi "error\|failed\|segmentation fault" "$temp_abundance_output"; then
                        if [ $abundance_exit_code -eq 0 ]; then
                            abundance_exit_code=1
                        fi
                    fi
                    
                    # Find the output file
                    local found_output=""
                    local possible_outputs=(
                        "$OUTPUT_DIR_ABS/${input_base}_abundance.csv"
                        "$OUTPUT_DIR_ABS/${input_base}.csv_abundance.csv"
                        "$OUTPUT_DIR_ABS/abundance.csv"
                        "$(dirname "$input_csv")/${input_base}_abundance.csv"
                        "$(dirname "$input_csv")/${input_base}.csv_abundance.csv"
                        "$CLARK_DIR/${input_base}_abundance.csv"
                    )
                    
                    for output_file in "${possible_outputs[@]}"; do
                        if [ -f "$output_file" ]; then
                            if [ "$(dirname "$output_file")" != "$OUTPUT_DIR_ABS" ]; then
                                cp "$output_file" "$abundance_output_file"
                            else
                                abundance_output_file="$output_file"
                            fi
                            found_output="$abundance_output_file"
                            break
                        fi
                    done
                    
                    # If extracting CSV data from the log
                    if [ -z "$found_output" ] && [ $abundance_exit_code -eq 0 ]; then
                        echo "Extracting abundance data from output log..."
                        local csv_start_line=$(grep -n "^Name,TaxID" "$temp_abundance_output" | head -1 | cut -d: -f1)
                        if [ -n "$csv_start_line" ]; then
                            sed -n "${csv_start_line},\$p" "$temp_abundance_output" > "$abundance_output_file"
                            found_output="$abundance_output_file"
                            echo "Extracted abundance data to: $abundance_output_file"
                        fi
                    fi
                    
                    rm -f "$temp_abundance_output"
                    
                    if [ -n "$found_output" ] && [ -f "$found_output" ]; then
                        ABUNDANCE_EXIT_CODE=$abundance_exit_code
                        ABUNDANCE_OUTPUT_FILE="$found_output"
                        return $abundance_exit_code
                    fi
                fi
                
                # Method 2: compute the abundance manually from the classification result CSV
                echo "Using manual abundance calculation from classification results..."
                echo "Start time: $(date)"
                echo "----------------------------------------"
                
                # Check the CSV file format
                if [ ! -f "$input_csv" ]; then
                    echo "ERROR: Classification file not found: $input_csv"
                    ABUNDANCE_EXIT_CODE=1
                    ABUNDANCE_OUTPUT_FILE=""
                    return 1
                fi
                
                # Use awk or Python to compute the abundance
                # First check the column structure of the CSV file
                local header=$(head -1 "$input_csv")
                echo "CSV header: $header"
                
                # Try to compute the abundance with Python (if available)
                if command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
                    local python_cmd=$(command -v python3 2>/dev/null || command -v python 2>/dev/null)
                    echo "Using Python for abundance calculation..."
                    
                    # Create a temporary Python script to compute the abundance
                    local py_script=$(mktemp)
                    cat > "$py_script" << 'PYEOF'
import sys
import csv
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the classification results
counts = defaultdict(int)
total_count = 0
header = None
data_rows = []

with open(input_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    
    # Find the index of the Count column
    count_col = None
    name_col = None
    for i, col in enumerate(header):
        if 'count' in col.lower() or col.lower() == 'count':
            count_col = i
        if 'name' in col.lower() or col.lower() == 'name':
            name_col = i
    
    if count_col is None or name_col is None:
        print("ERROR: Cannot find Count or Name column in CSV file", file=sys.stderr)
        sys.exit(1)
    
    # Read the data
    for row in reader:
        if len(row) > max(count_col, name_col):
            try:
                count = int(float(row[count_col]))
                name = row[name_col]
                counts[name] += count
                total_count += count
                data_rows.append(row)
            except (ValueError, IndexError):
                continue

# Compute the abundance and write the output file
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    
    # Write the header (adding the abundance column)
    new_header = list(header)
    if 'Proportion_All(%)' not in header:
        new_header.append('Proportion_All(%)')
    if 'Proportion_Classified(%)' not in header:
        new_header.append('Proportion_Classified(%)')
    writer.writerow(new_header)
    
    # Write the data
    for row in data_rows:
        if len(row) > max(count_col, name_col):
            try:
                count = int(float(row[count_col]))
                name = row[name_col]
                
                # Compute the abundance percentage
                prop_all = (count / total_count * 100) if total_count > 0 else 0
                prop_classified = (count / sum(counts.values()) * 100) if sum(counts.values()) > 0 else 0
                
                new_row = list(row)
                if 'Proportion_All(%)' not in header:
                    new_row.append(f"{prop_all:.6f}")
                if 'Proportion_Classified(%)' not in header:
                    new_row.append(f"{prop_classified:.6f}")
                
                writer.writerow(new_row)
            except (ValueError, IndexError):
                continue

print(f"Total sequences: {total_count}")
print(f"Total classified: {sum(counts.values())}")
PYEOF
                    
                    if $python_cmd "$py_script" "$input_csv" "$abundance_output_file" 2>&1; then
                        echo "----------------------------------------"
                        echo "End time: $(date)"
                        echo "Abundance calculation completed successfully!"
                        ABUNDANCE_EXIT_CODE=0
                        ABUNDANCE_OUTPUT_FILE="$abundance_output_file"
                        rm -f "$py_script"
                        return 0
                    else
                        echo "ERROR: Python abundance calculation failed"
                        rm -f "$py_script"
                    fi
                fi
                
                # If Python is not available, give an error message
                echo "ERROR: Python is not available for abundance calculation."
                echo "Please install Python 3 or ensure estimate_abundance.sh is available."
                echo "Classification results are still available at: $input_csv"
                ABUNDANCE_EXIT_CODE=1
                ABUNDANCE_OUTPUT_FILE=""
                return 1
            }
            
            # Run the abundance estimation
            if run_abundance_estimation "$CLASSIFICATION_CSV"; then
                echo "=========================================="
                echo "Abundance estimation completed successfully!"
                if [ -n "$ABUNDANCE_OUTPUT_FILE" ] && [ -f "$ABUNDANCE_OUTPUT_FILE" ]; then
                    echo "Abundance result file: $ABUNDANCE_OUTPUT_FILE"
                    echo ""
                    echo "Abundance file statistics:"
                    wc -l "$ABUNDANCE_OUTPUT_FILE"
                    echo ""
                    echo "Top 10 most abundant taxa:"
                    head -11 "$ABUNDANCE_OUTPUT_FILE" | column -t -s, 2>/dev/null || head -11 "$ABUNDANCE_OUTPUT_FILE"
                else
                    echo "WARNING: Abundance output file not found, but estimation completed."
                    echo "Check the output above for abundance data."
                fi
                echo "=========================================="
            else
                echo "=========================================="
                echo "WARNING: Abundance estimation failed or encountered errors."
                echo "Exit code: $ABUNDANCE_EXIT_CODE"
                echo "Classification results are still available at: $CLASSIFICATION_CSV"
                echo "=========================================="
            fi
        fi
    else
        echo ""
        echo "Note: Abundance estimation is disabled (ENABLE_ABUNDANCE=false)"
        echo "To enable, set ENABLE_ABUNDANCE=true before running the script"
    fi
else
    echo "=========================================="
    echo "Error: Virus classification failed! Exit code: $CLASSIFY_EXIT_CODE"
    
    # Show system resource information for diagnosis
    echo ""
    echo "System resource information:"
    if command -v free >/dev/null 2>&1; then
        echo "Memory usage:"
        free -h | head -2
    fi
    if command -v df >/dev/null 2>&1; then
        echo "Disk space in $TMPDIR:"
        df -h "$TMPDIR" | tail -1
    fi
    echo "Final thread count used: $CURRENT_THREADS"
    echo "Execution mode: $MODE"
    
    if [ "$ERROR_TYPE" = "segfault" ] || [ $CLASSIFY_EXIT_CODE -eq 139 ] || [ $CLASSIFY_EXIT_CODE -eq 11 ]; then
        echo ""
        echo "Segmentation fault occurred. Possible causes and solutions:"
        echo ""
        echo "1. Insufficient memory:"
        echo "   - Try reducing threads: export FORCE_THREADS=2 before running"
        echo "   - Request more memory in SLURM: --mem=256G"
        echo "   - Use express mode (mode 2) which uses less memory: MODE=2"
        echo ""
        echo "2. Database issues:"
        echo "   - Verify database integrity: check $CLARK_DB_DIR/Viruses"
        echo "   - Rebuild database if corrupted"
        echo "   - Ensure database matches k-mer size ($K_MER)"
        echo ""
        echo "3. Input file issues:"
        echo "   - Verify input files are not corrupted"
        echo "   - Check file sizes (very large files may cause issues)"
        echo "   - Try with smaller subset of data first"
        echo ""
        echo "4. CLARK binary compatibility:"
        echo "   - Check CLARK version compatibility"
        echo "   - Verify binary was compiled for this system architecture"
        echo ""
        echo "Recommended next steps:"
        echo "  - Retry with: FORCE_THREADS=2 MODE=2 sbatch $0"
        echo "  - Or manually reduce threads in the script (line ~150)"
    elif [ "$ERROR_TYPE" = "disk_space" ]; then
        echo ""
        echo "Disk space error occurred. Solutions:"
        echo "  - Set TMPDIR to a directory with more space:"
        echo "    export TMPDIR=/scratch/your_username/tmp"
        echo "  - Clean up temporary files:"
        echo "    rm -rf $TMPDIR/CLARKGZP*"
        echo "  - Use uncompressed input files (.fastq instead of .fastq.gz)"
        echo "  - Request more disk space or use a different partition"
    elif [ "$ERROR_TYPE" = "general" ]; then
        echo ""
        echo "General error occurred. Check CLARK output above for details."
        echo "Common issues:"
        echo "  - Incorrect file paths or permissions"
        echo "  - Database path issues"
        echo "  - CLARK configuration problems"
    fi
    echo "=========================================="
    exit $CLASSIFY_EXIT_CODE
fi

