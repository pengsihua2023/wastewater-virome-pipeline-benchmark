#!/bin/bash
#SBATCH --job-name=short_mag_708_2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=120:00:00
#SBATCH --output=short_mag_708_2_%j.out
#SBATCH --error=short_mag_708_2_%j.err

# ==============================================================================
# Illumina SHORT-READ metagenomics with the REAL nf-core/mag (Nextflow) pipeline.
# nf-core/mag natively supports short reads (assembly = MEGAHIT + metaSPAdes).
# Pipeline: fastp QC -> PhiX removal -> assembly -> Kraken2 taxonomy.
#
# Files used (all in the submit directory):
#   samplesheet_short.csv        input (sample,group,short_reads_1,short_reads_2,long_reads)
#   nfcore_mag_short.config      resources + Kraken2 Standard DB + skip binning/QC/GTDB-Tk
#
# Submit:   sbatch run_short_reads_mag_nfcore.sh
# Edited: 37, 39, 40
# ==============================================================================

set -Eeuo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1

# ---- CRLF guard ----
for f in run_short_reads_mag_nfcore.sh nfcore_mag_short.config samplesheet_short.csv; do
    if [ -f "$f" ] && grep -q $'\r' "$f"; then
        echo "ERROR: Windows CRLF detected in $f. Fix: dos2unix $f"; exit 1
    fi
done

# ------------------------------------------------------------------
# Config
# ------------------------------------------------------------------
SAMPLESHEET="${SAMPLESHEET:-samplesheet_short.csv}"
CONFIG="${CONFIG:-nfcore_mag_short.config}"
OUTDIR="${OUTDIR:-results_short_708_2}"
WORKDIR="${WORKDIR:-work_short_708_2}"
KRAKEN2_DB="${KRAKEN2_DB:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/k2_standard_20240605}"

echo "=========================================="
echo "nf-core/mag - Illumina short-read analysis"
echo "Samplesheet : $SAMPLESHEET"
echo "Config      : $CONFIG"
echo "Output      : $OUTDIR"
echo "Kraken2 DB  : $KRAKEN2_DB"
echo "Start       : $(date)"
echo "Job ID      : $SLURM_JOB_ID"
echo "=========================================="

[ -f "$SAMPLESHEET" ] || { echo "ERROR: samplesheet not found: $SAMPLESHEET"; exit 1; }
[ -f "$CONFIG" ]      || { echo "ERROR: config not found: $CONFIG"; exit 1; }
[ -d "$KRAKEN2_DB" ]  || { echo "ERROR: Kraken2 DB dir not found: $KRAKEN2_DB"; exit 1; }

# ------------------------------------------------------------------
# Conda / Nextflow environment
# ------------------------------------------------------------------
module load Miniforge3/24.11.3-0 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"

CONDA_ENV_PATH="/home/$USER/.conda/envs/nextflow_env"
if [ -d "$CONDA_ENV_PATH" ]; then
    conda activate "$CONDA_ENV_PATH" || export PATH="$CONDA_ENV_PATH/bin:$PATH"
    echo "Activated conda env: $CONDA_ENV_PATH"
else
    echo "ERROR: conda env not found at $CONDA_ENV_PATH"; exit 1
fi

# The login env sets PYTHONPATH/PYTHONHOME (user-py39) which breaks conda python3.
unset PYTHONPATH PYTHONHOME

# Java check (Nextflow needs 11-23)
JAVA_MAJOR=$(java -version 2>&1 | head -1 | grep -oE '"([0-9]+)' | grep -oE '[0-9]+' | head -1)
if [ -n "$JAVA_MAJOR" ] && { [ "$JAVA_MAJOR" -lt 11 ] || [ "$JAVA_MAJOR" -gt 23 ]; }; then
    echo "ERROR: Java $JAVA_MAJOR incompatible with Nextflow (needs 11-23)"; exit 1
fi
echo "Java: ${JAVA_MAJOR:-unknown}"

# ------------------------------------------------------------------
# Nextflow environment variables
# ------------------------------------------------------------------
export NXF_HOME="/scratch/$USER/Meta-genome-data-analysis/Apptainer/nextflow"
mkdir -p "$NXF_HOME"
NXF_USER_HOME="/scratch/$USER/Meta-genome-data-analysis/Apptainer/nextflow_user_home"
mkdir -p "$NXF_USER_HOME"
cat > "$NXF_USER_HOME/.gitconfig" << 'EOF_GIT'
[core]
    autocrlf = false
EOF_GIT
export NXF_OPTS="-Xms1g -Xmx4g -Dnf.validation.lenient=true -Duser.home=${NXF_USER_HOME}"
export APPTAINER_NO_MOUNT="/lscratch"
export SINGULARITY_NO_MOUNT="/lscratch"
export NXF_SINGULARITY_CACHEDIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/singularity"
export NXF_APPTAINER_CACHEDIR="$NXF_SINGULARITY_CACHEDIR"

# ------------------------------------------------------------------
# Pull nf-core/mag r3.1.0 (idempotent)
# ------------------------------------------------------------------
# A previous long-read attempt patched the cached schema_input.json, so the
# cached repo has "uncommitted changes" and `nextflow pull` refuses. Discard
# those local edits first (short reads don't need that patch), then pull.
ASSETS="${NXF_HOME}/assets/nf-core/mag"
if [ -d "$ASSETS/.git" ]; then
    echo "Resetting cached nf-core/mag (discarding local edits)..."
    git -C "$ASSETS" checkout -- . 2>/dev/null || true
    git -C "$ASSETS" clean -fdx 2>/dev/null || true
fi
echo "Pulling nf-core/mag r3.1.0 ..."
nextflow pull nf-core/mag -r 3.1.0 || echo "WARN: pull failed; using existing cached assets at $ASSETS"

# ------------------------------------------------------------------
# Kraken2 DB for mag's built-in READ-LEVEL classification.
# nf-core/mag 3.1.0's --kraken2_db expects a .tar.gz archive (it untars it in a
# KRAKEN2_DB_PREPARATION step); a plain directory or placeholder does NOT work
# (that was why the old placeholder + ext.args hack produced 0 Kraken2 tasks).
# To avoid copying/taring the ~140GB standard DB, build a TINY tarball that
# contains SYMLINKS to the real *.k2d files. kraken2 resolves them at run time
# because the config bind-mounts /scratch into the container.
# Fallback (if mag/kraken2 cannot follow the symlinks): make a real archive with
#   tar -czf "$KRAKEN2_DB_TGZ" -C "$KRAKEN2_DB" .
KRAKEN2_DB_TGZ="${NXF_HOME}/kraken2_db_links.tar.gz"
if [ ! -f "$KRAKEN2_DB_TGZ" ]; then
    echo "Building lightweight Kraken2 DB tarball (symlinks) for mag..."
    LINKROOT="${NXF_HOME}/kraken2_db_links"
    rm -rf "$LINKROOT"
    mkdir -p "$LINKROOT/database"
    _linked=0
    for f in hash.k2d opts.k2d taxo.k2d seqid2taxid.map; do
        if [ -e "$KRAKEN2_DB/$f" ]; then
            ln -sf "$KRAKEN2_DB/$f" "$LINKROOT/database/$f"
            [ "$f" != "seqid2taxid.map" ] && _linked=$((_linked+1))
        fi
    done
    if [ "$_linked" -lt 3 ]; then
        echo "ERROR: $KRAKEN2_DB is missing hash.k2d/opts.k2d/taxo.k2d — not a valid Kraken2 DB."
        exit 1
    fi
    # -C ... database  ->  archive stores the SYMLINKS themselves (no -h), tiny.
    tar -czf "$KRAKEN2_DB_TGZ" -C "$LINKROOT" database
    echo "  Created: $KRAKEN2_DB_TGZ ($(du -h "$KRAKEN2_DB_TGZ" | cut -f1))"
fi

# Remove stale report files so -resume does not error on overwrite
rm -f "$OUTDIR/execution_report.html" "$OUTDIR/execution_timeline.html" 2>/dev/null

# ------------------------------------------------------------------
# nf-core/mag command
# ------------------------------------------------------------------
NF_CMD=(
    nextflow run nf-core/mag
    -r 3.1.0
    -profile apptainer
    -c "$CONFIG"
    --validationSkipDuplicateCheck
    --input "$SAMPLESHEET"
    --outdir "$OUTDIR"
    --max_cpus 32
    --max_memory '256.GB'
    --max_time '120.h'
    --coassemble_group
    --megahit_options '--min-contig-len 1000'
    --spades_options '--meta --only-assembler'
    --fastp_save_trimmed_fail
    --fastp_qualified_quality 20
    --save_clipped_reads
    --kraken2_db "$KRAKEN2_DB_TGZ"
    --skip_binning
    --skip_binqc
    --skip_gtdbtk
    --skip_prokka
    --skip_quast
    -work-dir "$WORKDIR"
    # -resume           # Commented out for fresh run; uncomment to resume from previous work-dir
    -with-report "$OUTDIR/execution_report.html"
    -with-timeline "$OUTDIR/execution_timeline.html"
)

echo "Command:"
printf '%q ' "${NF_CMD[@]}"; echo ""
echo "=========================================="

"${NF_CMD[@]}"
EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "nf-core/mag completed successfully!"
    echo "Results : $OUTDIR/"

    # Only report directories that actually exist
    [ -d "$OUTDIR/Assembly" ] && echo "  Assemblies : $OUTDIR/Assembly/" || echo "  Assemblies : (not generated)"
    [ -d "$OUTDIR/Taxonomy/kraken2" ] && echo "  Kraken2    : $OUTDIR/Taxonomy/kraken2/" || echo "  Kraken2    : (not generated - check if --skip_kraken2 was set or Kraken2 is enabled by default)"
    [ -d "$OUTDIR/MultiQC" ] && echo "  MultiQC    : $OUTDIR/MultiQC/" || echo "  MultiQC    : (not generated)"
else
    echo "ERROR: nf-core/mag failed with exit code $EXIT_CODE"
fi
echo "End: $(date)"
echo "=========================================="
exit $EXIT_CODE
