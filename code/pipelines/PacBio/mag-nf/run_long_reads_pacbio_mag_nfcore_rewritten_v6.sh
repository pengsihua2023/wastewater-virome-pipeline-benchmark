#!/usr/bin/env bash
#SBATCH --job-name=mag_pacbio_sam0
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=120:00:00
#SBATCH --output=mag_pacbio_sam0_%j.out
#SBATCH --error=mag_pacbio_sam0_%j.err

# ==============================================================================
# nf-core/mag PacBio HiFi long-read-only metagenome assembly workflow
# v6 fixes:
#   1) uses -params-file JSON so booleans/numbers stay typed correctly
#   2) avoids --skip_binqc + default --run_busco=true conflict by setting run_busco=false
#   3) keeps only infrastructure/resource settings in -c config
# Edited:28(samplesheet),out:31;work:32
# ==============================================================================

set -Eeuo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    cd "$SLURM_SUBMIT_DIR"
fi

PIPELINE_VERSION="${PIPELINE_VERSION:-5.4.2}"
NEXTFLOW_VERSION="${NEXTFLOW_VERSION:-26.04.4}"
SAMPLESHEET="${SAMPLESHEET:-samplesheet_pacbio_rewritten.csv}"
CONFIG="${CONFIG:-nfcore_mag_pacbio_rewritten_v6.config}"
PARAMS_FILE="${PARAMS_FILE:-params_pacbio_mag_v6.json}"
OUTDIR="${OUTDIR:-results_pacbio_hifi_sample_0}"
WORKDIR="${WORKDIR:-work_pacbio_hifi_sample_0}"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-/home/$USER/.conda/envs/nextflow_env}"
NXF_BASE="${NXF_BASE:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/nextflow}"
CONTAINER_CACHE="${CONTAINER_CACHE:-/scratch/$USER/Meta-genome-data-analysis/Apptainer/singularity}"

COASSEMBLE_GROUP="${COASSEMBLE_GROUP:-false}"
# metaMDBG is built for large real HiFi metagenomes (millions of reads) and
# crashes at the "Converting reads to minimizers" stage on the tiny simulated
# CAMISIM sample (~72 reads), aborting the whole run. Default to OFF and let
# Flye (metaFlye) do the long-read assembly — the same assembler the sourmash
# PacBio pipeline uses on this data. Set RUN_METAMDBG=true to re-enable.
RUN_METAMDBG="${RUN_METAMDBG:-false}"
RUN_VIRUS_IDENTIFICATION="${RUN_VIRUS_IDENTIFICATION:-false}"
GENOMAD_DB="${GENOMAD_DB:-}"

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { echo "ERROR: $*" >&2; exit 1; }

version_ge() {
    python3 - "$1" "$2" <<'PYVER'
import re, sys
def nums(v):
    parts = re.findall(r'\d+', v)[:3]
    return tuple(int(x) for x in (parts + ['0','0','0'])[:3])
sys.exit(0 if nums(sys.argv[1]) >= nums(sys.argv[2]) else 1)
PYVER
}

get_nextflow_version() {
    local nf="$1"
    [[ -x "$nf" ]] || return 1
    "$nf" -version 2>&1 | awk '/version/ {for(i=1;i<=NF;i++) if($i ~ /^[0-9]+\.[0-9]+/) {print $i; exit}}'
}

install_private_nextflow() {
    local nf_bin="$1"
    local tmp_dir="${NXF_BASE}/tmp/install_nextflow_${NEXTFLOW_VERSION}_$$"
    mkdir -p "$(dirname "$nf_bin")" "$tmp_dir"
    log "Installing private Nextflow ${NEXTFLOW_VERSION} to ${nf_bin}"
    (
        set -Eeuo pipefail
        cd "$tmp_dir"
        if command -v curl >/dev/null 2>&1; then
            curl -fsSL https://get.nextflow.io -o get.nextflow.io
        elif command -v wget >/dev/null 2>&1; then
            wget -qO get.nextflow.io https://get.nextflow.io
        else
            echo "ERROR: curl/wget not found. Install Nextflow manually, or load a module with Nextflow >=25.04.0." >&2
            exit 1
        fi
        chmod +x get.nextflow.io
        NXF_VER="$NEXTFLOW_VERSION" NXF_HOME="$NXF_BASE" bash get.nextflow.io >/dev/null
        [[ -f nextflow ]] || { echo "ERROR: Nextflow installer did not create $tmp_dir/nextflow" >&2; exit 1; }
        chmod +x nextflow
        mv -f nextflow "$nf_bin"
    )
    rm -rf "$tmp_dir"
    [[ -x "$nf_bin" ]] || die "Nextflow executable was not created: $nf_bin"
}

select_nextflow() {
    local min_version="25.04.0"
    local nf_bin="${NXF_BASE}/bin/nextflow"
    mkdir -p "${NXF_BASE}/bin"

    local private_ver=""
    private_ver="$(get_nextflow_version "$nf_bin" || true)"
    if [[ -n "$private_ver" ]] && version_ge "$private_ver" "$min_version"; then
        echo "$nf_bin"
        return 0
    fi

    local system_nf=""
    system_nf="$(command -v nextflow || true)"
    if [[ -n "$system_nf" ]]; then
        local system_ver=""
        system_ver="$(get_nextflow_version "$system_nf" || true)"
        if [[ -n "$system_ver" ]] && version_ge "$system_ver" "$min_version"; then
            echo "$system_nf"
            return 0
        fi
        log "Found old Nextflow at ${system_nf}: ${system_ver:-unknown}; nf-core/mag ${PIPELINE_VERSION} needs >=${min_version}."
    fi

    install_private_nextflow "$nf_bin"
    private_ver="$(get_nextflow_version "$nf_bin" || true)"
    [[ -n "$private_ver" ]] || die "Installed Nextflow but could not read its version: $nf_bin"
    version_ge "$private_ver" "$min_version" || die "Installed Nextflow $private_ver is too old; need >=$min_version"
    echo "$nf_bin"
}

log "nf-core/mag PacBio HiFi long-read-only run"
log "Pipeline version : ${PIPELINE_VERSION}"
log "Nextflow target  : ${NEXTFLOW_VERSION}"
log "Samplesheet      : ${SAMPLESHEET}"
log "Config           : ${CONFIG}"
log "Params file      : ${PARAMS_FILE}"
log "Outdir           : ${OUTDIR}"
log "Workdir          : ${WORKDIR}"
log "Job ID           : ${SLURM_JOB_ID:-not_slurm}"

[[ -f "$SAMPLESHEET" ]] || die "Samplesheet not found: $SAMPLESHEET"
[[ -f "$CONFIG" ]] || die "Config not found: $CONFIG"

for f in "$0" "$SAMPLESHEET" "$CONFIG"; do
    if [[ -f "$f" ]] && grep -q $'\r' "$f"; then
        die "Windows CRLF detected in $f. Run: dos2unix $f"
    fi
done

header="$(head -n 1 "$SAMPLESHEET")"
expected_header="sample,run,group,long_reads,long_reads_platform"
[[ "$header" == "$expected_header" ]] || die "Bad samplesheet header. Expected: $expected_header ; Got: $header"
awk -F',' 'NR>1 { if (NF != 5) { printf("Bad CSV column count on line %d: expected 5, got %d\n", NR, NF); exit 1 } }' "$SAMPLESHEET"
awk -F',' 'NR>1 { if ($5 != "PACBIO_HIFI" && $5 != "PACBIO_CLR" && $5 != "OXFORD_NANOPORE" && $5 != "OXFORD_NANOPORE_HQ") { printf("Invalid long_reads_platform on line %d: %s\n", NR, $5); exit 1 } }' "$SAMPLESHEET"
while IFS=',' read -r sample run group long_reads platform; do
    [[ "$sample" == "sample" ]] && continue
    [[ -n "$sample" && -n "$run" && -n "$group" && -n "$long_reads" && -n "$platform" ]] || die "Empty field in samplesheet for sample: ${sample:-unknown}"
    [[ "$long_reads" == *.fastq.gz || "$long_reads" == *.fq.gz ]] || die "Long reads must be compressed FASTQ (.fastq.gz or .fq.gz): $long_reads"
    [[ -f "$long_reads" ]] || die "Long-read FASTQ not found: $long_reads"
done < "$SAMPLESHEET"

# ------------------------------------------------------------------
# Deduplicate long-read IDs.
# CAMISIM-simulated HiFi reads reuse read names (e.g. "S1/1/ccs"), which makes
# mag's FILTLONG abort with: "Error: duplicate read name: S1/1/ccs". Rewrite
# every read header to a unique ID and repoint the samplesheet at the cleaned
# FASTQs. Idempotent: an existing dedup file is reused.
# ------------------------------------------------------------------
DEDUP_DIR="${DEDUP_DIR:-dedup_reads}"
DEDUP_SAMPLESHEET="${DEDUP_SAMPLESHEET:-samplesheet_pacbio_dedup.csv}"
mkdir -p "$DEDUP_DIR"
head -n 1 "$SAMPLESHEET" > "$DEDUP_SAMPLESHEET"
while IFS=',' read -r sample run group long_reads platform; do
    [[ "$sample" == "sample" || -z "$sample" ]] && continue
    base="$(basename "$long_reads")"; base="${base%.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    dedup_out="${DEDUP_DIR}/${base}.dedup.fastq.gz"
    if [[ -s "$dedup_out" ]]; then
        log "Dedup file exists, reusing: $dedup_out"
    else
        log "Deduplicating read IDs: $long_reads -> $dedup_out"
        # Rewrite only the FASTQ header line (every 4th line) to a unique name.
        zcat "$long_reads" | awk 'BEGIN{c=0} NR%4==1{c++; print "@read_"c; next} {print}' | gzip > "$dedup_out"
    fi
    printf '%s,%s,%s,%s,%s\n' "$sample" "$run" "$group" "$(readlink -f "$dedup_out")" "$platform" >> "$DEDUP_SAMPLESHEET"
done < "$SAMPLESHEET"
# Export so the params-JSON generator (reads os.environ) uses the cleaned sheet.
export SAMPLESHEET="$DEDUP_SAMPLESHEET"
log "Using deduplicated samplesheet: $SAMPLESHEET"

module load Miniforge3/24.11.3-0 2>/dev/null || true
if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
if [[ -d "$CONDA_ENV_PATH" ]]; then
    conda activate "$CONDA_ENV_PATH" || export PATH="$CONDA_ENV_PATH/bin:$PATH"
else
    log "WARN: Conda env not found: $CONDA_ENV_PATH. Continuing with current environment."
fi
unset PYTHONPATH PYTHONHOME
command -v java >/dev/null 2>&1 || die "java not found. Load Java 11-23 or activate a conda env that provides Java."
command -v python3 >/dev/null 2>&1 || die "python3 not found; needed for version checks."
JAVA_MAJOR="$(java -version 2>&1 | awk -F'\"' '/version/ {print $2}' | awk -F. '{print ($1 == "1") ? $2 : $1}')"
if [[ -n "$JAVA_MAJOR" ]] && (( JAVA_MAJOR < 11 || JAVA_MAJOR > 23 )); then
    die "Java $JAVA_MAJOR is incompatible with Nextflow. Use Java 11-23."
fi

export NXF_HOME="$NXF_BASE"
mkdir -p "$NXF_HOME" "$CONTAINER_CACHE" "$OUTDIR" "$WORKDIR"
NXF_USER_HOME="${NXF_BASE}/nextflow_user_home"
mkdir -p "$NXF_USER_HOME"
cat > "${NXF_USER_HOME}/.gitconfig" <<'GITCONFIG'
[core]
    autocrlf = false
GITCONFIG

export NXF_OPTS="-Xms1g -Xmx4g -Duser.home=${NXF_USER_HOME}"
# Keep parser v1 for nf-core/mag 5.4.2 under Nextflow 26.x.
export NXF_SYNTAX_PARSER="v1"
export APPTAINER_NO_MOUNT="/lscratch"
export SINGULARITY_NO_MOUNT="/lscratch"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE"

NEXTFLOW_BIN="$(select_nextflow)"
log "Using Nextflow executable: $NEXTFLOW_BIN"
log "Nextflow version: $($NEXTFLOW_BIN -version 2>&1 | head -n 2 | tr '\n' ' ')"
log "Java major version: ${JAVA_MAJOR:-unknown}"
log "NXF_SYNTAX_PARSER: ${NXF_SYNTAX_PARSER}"

# Write params as JSON to preserve boolean and numeric types for nf-schema.
python3 - <<PYJSON
import json, os
params = {
    "input": os.environ.get("SAMPLESHEET", "${SAMPLESHEET}"),
    "outdir": os.environ.get("OUTDIR", "${OUTDIR}"),
    "max_cpus": 32,
    "max_memory": "256.GB",
    "max_time": "120.h",
    "skip_shortread_qc": True,
    "skip_binning": True,
    "skip_binqc": True,
    "run_busco": False,
    "skip_gtdbtk": True,
    "skip_prokka": True,
    "skip_metaeuk": True,
    "skip_quast": True,
    "skip_ale": True,
    "save_filtered_longreads": True,
    "coassemble_group": os.environ.get("COASSEMBLE_GROUP", "${COASSEMBLE_GROUP}").lower() == "true",
}
if os.environ.get("RUN_METAMDBG", "${RUN_METAMDBG}").lower() != "true":
    params["skip_metamdbg"] = True
if os.environ.get("RUN_VIRUS_IDENTIFICATION", "${RUN_VIRUS_IDENTIFICATION}").lower() == "true":
    genomad_db = os.environ.get("GENOMAD_DB", "${GENOMAD_DB}")
    if not genomad_db:
        raise SystemExit("RUN_VIRUS_IDENTIFICATION=true requires GENOMAD_DB")
    params["run_virus_identification"] = True
    params["genomad_db"] = genomad_db
with open("${PARAMS_FILE}", "w") as f:
    json.dump(params, f, indent=2)
    f.write("\n")
PYJSON

log "Generated typed params file: ${PARAMS_FILE}"
cat "${PARAMS_FILE}" >&2

log "Pulling nf-core/mag ${PIPELINE_VERSION}"
"$NEXTFLOW_BIN" pull nf-core/mag -r "$PIPELINE_VERSION"
rm -f "$OUTDIR/execution_report.html" "$OUTDIR/execution_timeline.html" "$OUTDIR/execution_trace.txt" 2>/dev/null || true

NF_CMD=(
    "$NEXTFLOW_BIN" run nf-core/mag
    -r "$PIPELINE_VERSION"
    -profile apptainer
    -c "$CONFIG"
    -params-file "$PARAMS_FILE"
    -work-dir "$WORKDIR"
    -with-report "$OUTDIR/execution_report.html"
    -with-timeline "$OUTDIR/execution_timeline.html"
    -with-trace "$OUTDIR/execution_trace.txt"
)

log "Command:"
printf '%q ' "${NF_CMD[@]}" >&2
echo >&2

"${NF_CMD[@]}"
exit_code=$?
if [[ $exit_code -eq 0 ]]; then
    log "nf-core/mag completed successfully. Results: $OUTDIR"
else
    log "nf-core/mag failed with exit code $exit_code"
fi
exit "$exit_code"
