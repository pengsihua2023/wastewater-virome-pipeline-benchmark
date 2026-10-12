#!/bin/bash
#SBATCH --job-name=CLARK_sam2_pacbio
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          # raised from 4 -> faster classify; lower if your partition limits it
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=CLARK_sam2_pacbio_%j.out
#SBATCH --error=CLARK_sam2_pacbio_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ------------------------------------------
# Configuration   Edited lines 16, 17, 18
# ------------------------------------------
INPUT_FILE="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/pacbio-708-genomes/pacbio_hifi_reads_sample_2.fastq.gz"
OUTPUT_PREFIX="pacbio_hifi_708_2_virus"
OUTPUT_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/CLARK/Pacbio-708/results_long-pacbio-sample-2"
CLARK_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/CLARK"

# CLARK custom-DB directory (where set_targets.sh built the targets index).
# This is the SAME DIR_DB used in set_targets.sh, confirmed from $CLARK_DIR/.settings:
#   -T <DIR_DB>/targets.txt
#   -D <DIR_DB>/viruses_0/      (viruses-only, species level -> RANK 0)
DIR_DB="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/CLARK_DB"

K_MER=20     # MUST match the k used to BUILD the database (viruses_0 was built at k=20).
             # Classifying at a different k than the DB was built with yields ALL-NA results.
             # To use a different k (e.g. 31 for HiFi specificity), you must REBUILD the DB
             # at that k (SETUP_TARGETS=1) AND remove the old viruses_0 index first.
MODE=0       # KEEP 0 (full): only full mode emits gamma / confidence / 2nd_assignment needed for downstream filtering. Do not change.

# Gamma is CLARK's per-read match strength = (score1+score2)/#kmers-in-read, range 0-1.
# With a virus-only DB on a mostly-bacterial sample, bacterial reads get force-assigned
# to giant viruses but with VERY low gamma (<0.05); real viral reads have gamma ~0.26+.
# Filtering at GAMMA_MIN cleanly removes that bacterial noise: precision 2% -> ~90%,
# sensitivity stays ~88% (validated on pacbio sample_0 vs species_list_final.tsv).
#   0.2  = recommended (precision 90% / sens 88%, keeps the dominant astrovirus)
#   0.3  = higher precision (91% / 85%)
#   DO NOT exceed 0.4 -> starts dropping real low-gamma viruses incl. the dominant species.
GAMMA_MIN=0.2

# Set to 1 to (re)build the VIRUS targets database (species level) before classifying.
# Normally keep 0: the existing DIR_DB/viruses_0 index is reused as-is.
# Only set to 1 if you change k-mer/rank or the viral reference set in DIR_DB.
SETUP_TARGETS=0

echo "=========================================="
echo "  CLARK Virus Classification (PacBio HiFi)"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

echo "Configuration:"
echo "  CLARK directory:  $CLARK_DIR"
echo "  DB directory:     $DIR_DB"
echo "  Input file:       $INPUT_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Output prefix:    $OUTPUT_PREFIX"
echo "  k-mer size:       $K_MER"
echo "  Mode:             $MODE"
echo "  Rebuild targets:  $SETUP_TARGETS"
echo ""

# ------------------------------------------
# 2. Validate CLARK installation
# ------------------------------------------
CLASSIFY_SCRIPT="$CLARK_DIR/classify_metagenome.sh"
SET_TARGETS_SCRIPT="$CLARK_DIR/set_targets.sh"
ESTIMATE_SCRIPT="$CLARK_DIR/estimate_abundance.sh"

if [ ! -f "$CLASSIFY_SCRIPT" ]; then
    echo "ERROR: classify_metagenome.sh not found: $CLASSIFY_SCRIPT"
    exit 1
fi
echo "Classify script: $CLASSIFY_SCRIPT"

# ------------------------------------------
# 3. Validate input file
# ------------------------------------------
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi
echo "Input file size: $(du -h "$INPUT_FILE" | cut -f1)"
echo ""

# ------------------------------------------
# 4. Create output directory
# ------------------------------------------
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------
# 4b. (Optional) Build VIRUS + BACTERIA decoy targets at species level
# ------------------------------------------
if [ "$SETUP_TARGETS" -eq 1 ]; then
    echo "Setting CLARK targets (viruses, species level)..."
    # Order/args: set_targets.sh <DIR_DB> <db1> [db2 ...] --species
    # Virus-only -> builds DIR_DB/viruses_0, matching the existing .settings.
    "$SET_TARGETS_SCRIPT" "$DIR_DB" viruses --species
    if [ $? -ne 0 ]; then
        echo "ERROR: set_targets.sh failed"; exit 1
    fi
    echo "Targets set."
    echo ""
fi

# ------------------------------------------
# 5. Decompress .gz if needed (CLARK does not read .gz)
# ------------------------------------------
if [[ "$INPUT_FILE" == *.gz ]]; then
    DECOMPRESSED="${OUTPUT_DIR}/$(basename "$INPUT_FILE" .gz)"
    echo "Decompressing input file..."
    gunzip -c "$INPUT_FILE" > "$DECOMPRESSED"
    echo "Decompressed to: $DECOMPRESSED ($(du -h "$DECOMPRESSED" | cut -f1))"
    INPUT_FILE="$DECOMPRESSED"
    echo ""
fi

# ------------------------------------------
# 6. Run CLARK classification
# ------------------------------------------
OUTPUT_PATH="$OUTPUT_DIR/$OUTPUT_PREFIX"

echo "Running CLARK classification..."
echo "Command: $CLASSIFY_SCRIPT -O $INPUT_FILE -R $OUTPUT_PATH -m $MODE -n $SLURM_CPUS_PER_TASK -k $K_MER"
echo ""

"$CLASSIFY_SCRIPT" \
    -O "$INPUT_FILE" \
    -R "$OUTPUT_PATH" \
    -m "$MODE" \
    -n "$SLURM_CPUS_PER_TASK" \
    -k "$K_MER"

EXIT_CODE=$?

# ------------------------------------------
# 7. Abundance estimation -> clean table with TaxID + counts + proportion
# ------------------------------------------
RESULT_CSV="${OUTPUT_PATH}.csv"
if [ $EXIT_CODE -eq 0 ] && [ -f "$ESTIMATE_SCRIPT" ] && [ -f "$RESULT_CSV" ]; then
    echo ""
    # getAbundance reads the assignment ONLY from a 3-column (default-mode) file:
    #   Object_ID,Length,Assignment   (taxid in the LAST column)
    # Our classify output is full mode (8 cols), so we must extract cols 1,2,4
    # (Object_ID, Length, 1st_assignment) into a 3-column file first, else every
    # read -> UNKNOWN.  -D must come before -F (-F is variadic, swallows trailing args).

    # (a) UNFILTERED abundance: every classified read (precision is low here - bacterial
    #     reads pollute viral calls; kept for reference / full picture).
    echo "Estimating abundance (unfiltered, all classified reads)..."
    ABUND_ALL="${OUTPUT_PATH}.assignments.csv"
    awk -F, 'NR==1{print "Object_ID,Length,Assignment";next}{print $1","$2","$4}' \
        "$RESULT_CSV" > "$ABUND_ALL"
    "$ESTIMATE_SCRIPT" -D "$DIR_DB" -F "$ABUND_ALL" \
        > "${OUTPUT_PATH}.abundance.csv" 2> "${OUTPUT_PATH}.abundance.log"
    echo "Wrote: ${OUTPUT_PATH}.abundance.csv"

    # (b) GAMMA-FILTERED abundance: keep only reads with gamma (col 3) >= GAMMA_MIN
    #     and a real assignment. This is the high-precision viral detection table.
    echo "Estimating abundance (gamma >= ${GAMMA_MIN}, high-precision viral set)..."
    ABUND_G="${OUTPUT_PATH}.g${GAMMA_MIN}.assignments.csv"
    awk -F, -v G="$GAMMA_MIN" \
        'NR==1{print "Object_ID,Length,Assignment";next} ($3+0)>=G && $4!="NA"{print $1","$2","$4}' \
        "$RESULT_CSV" > "$ABUND_G"
    KEPT=$(( $(wc -l < "$ABUND_G") - 1 ))
    echo "  reads passing gamma>=${GAMMA_MIN}: $KEPT"
    "$ESTIMATE_SCRIPT" -D "$DIR_DB" -F "$ABUND_G" \
        > "${OUTPUT_PATH}.abundance.g${GAMMA_MIN}.csv" 2> "${OUTPUT_PATH}.abundance.g${GAMMA_MIN}.log"
    echo "Wrote: ${OUTPUT_PATH}.abundance.g${GAMMA_MIN}.csv  <-- use this for viral F1"
fi

echo ""
echo "=========================================="
echo "  Classification Results"
echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Classification completed successfully!"
    echo ""
    echo "Output files:"
    ls -lh "$OUTPUT_DIR/${OUTPUT_PREFIX}"* 2>/dev/null
    echo ""
    if [ -f "$RESULT_CSV" ]; then
        TOTAL_READS=$(( $(wc -l < "$RESULT_CSV") - 1 ))
        # FIX: assignment is column 4 (1st_assignment), not column 2 (Length).
        CLASSIFIED=$(awk -F',' 'NR>1 && $4 != "NA"' "$RESULT_CSV" | wc -l)
        echo "Summary:"
        echo "  Total reads:      $TOTAL_READS"
        echo "  Classified reads: $CLASSIFIED"
        echo ""
        echo "First 10 lines:"
        head -10 "$RESULT_CSV"
    fi
else
    echo "ERROR: Classification failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi
