#!/bin/bash
#SBATCH --job-name=GOTTCHA2-sam-2
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Pacbio_GOTTCHA2_sample_2_%j.out
#SBATCH --error=Pacbio_GOTTCHA2_sample_2_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "Metagenome Viral Classification Workflow (PacBio HiFi Mode)"

echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Node: ${SLURM_NODELIST:-NA}"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

echo "1. Setting up environment..."
module load Miniforge3/24.11.3-0

USER_CONDA_ENV="/home/$USER/.conda/envs/gottcha2_env"
source activate "$USER_CONDA_ENV"

# ============================================================
# GOTTCHA2 database settings
# ============================================================
# Your extracted GOTTCHA2 pan database contains:
#   gottcha_db.species.fna.mmi
#   gottcha_db.species.fna.stats
#   gottcha_db.species.fna.tax.tsv
#
# Therefore, -d should point to the database PREFIX without .mmi.
# GOTTCHA2 will automatically append .mmi internally.
# The real index file is ${DB_PREFIX}.mmi.
# Edited: 52, 53, 
# ============================================================

DB_DIR="/scratch/$USER/Meta-genome-data-analysis/Apptainer/databases/GOTTCHA2"

DB_PREFIX="${DB_DIR}/gottcha_db.species.fna"
DB_INDEX="${DB_PREFIX}.mmi"
DB_TAX="${DB_PREFIX}.tax.tsv"
DB_STATS="${DB_PREFIX}.stats"

INPUT="/scratch/$USER/Meta-genome-data-analysis/Apptainer/simulated-data-validation/data/pacbio-708-genomes/pacbio_hifi_reads_sample_2.fastq.gz"
OUTPUT_PREFIX="result_GOTTCHA2_pacbio_sample_2"

GOTTCHA2_BIN="/home/$USER/.conda/envs/gottcha2_env/bin/gottcha2.py"

echo "Checking input and database files..."
echo "Database directory: $DB_DIR"
echo "Database prefix:    $DB_PREFIX"
echo "Database index:     $DB_INDEX"
echo "Taxonomy file:      $DB_TAX"
echo "Stats file:         $DB_STATS"
echo "Input FASTQ:        $INPUT"
echo ""

ls -lh "$DB_INDEX" "$DB_TAX" "$DB_STATS" "$INPUT"

if [ ! -s "$DB_INDEX" ]; then
    echo "ERROR: GOTTCHA2 database index not found or empty: $DB_INDEX" >&2
    exit 1
fi

if [ ! -s "$DB_TAX" ]; then
    echo "ERROR: GOTTCHA2 taxonomy file not found or empty: $DB_TAX" >&2
    exit 1
fi

if [ ! -s "$DB_STATS" ]; then
    echo "ERROR: GOTTCHA2 database stats file not found or empty: $DB_STATS" >&2
    exit 1
fi

if [ ! -s "$INPUT" ]; then
    echo "ERROR: Input FASTQ file not found or empty: $INPUT" >&2
    exit 1
fi

if [ ! -x "$GOTTCHA2_BIN" ]; then
    echo "ERROR: GOTTCHA2 executable not found or not executable: $GOTTCHA2_BIN" >&2
    echo "Please check your conda environment or run: which gottcha2.py" >&2
    exit 1
fi

echo ""
echo "GOTTCHA2 executable:"
"$GOTTCHA2_BIN" --help | head -20 || true
echo ""

# Patch gottcha2.py to support PacBio HiFi data via -xm map-hifi.
# Four fixes are included:
#   1. 'map-hifi' not in -xm argparse choices  -> add it
#   2. AS tag parsed as raw string in remove_multiple_hits -> extract integer
#   3. parse() crashes on SAM lines with < 11 fields -> guard before temp[5]
#   4. worker() skips None returns from parse()
#
# These patches are written to be idempotent, so running the SLURM file
# multiple times should not repeatedly corrupt gottcha2.py.
echo "2. Patching gottcha2.py for map-hifi support..."
python3 - "$GOTTCHA2_BIN" << 'PYEOF'
import sys
import re

path = sys.argv[1]

with open(path) as f:
    src = f.read()

report = []
failed = []

# ── Patch 1: add map-hifi to -xm / --presetx argparse choices ─────────
def _add_map_hifi(m):
    s = m.group(0)
    if 'map-hifi' in s:
        return s
    q = "'" if "'map-ont'" in s else '"'
    return re.sub(
        rf"{q}map-ont{q}(\s*)\]",
        rf"{q}map-ont{q}\1, {q}map-hifi{q}]",
        s
    )

new_src = re.sub(r"choices\s*=\s*\[[^\]]*map-ont[^\]]*\]", _add_map_hifi, src)

if new_src != src:
    src = new_src
    report.append("Patch 1 applied: map-hifi added to -xm choices.")
elif 'map-hifi' in src and ('presetx' in src or '-xm' in src):
    report.append("Patch 1 already applied.")
else:
    failed.append("Patch 1 FAILED: choices list not found; map-hifi cannot be used.")

# ── Patch 2: fix AS tag value extraction in remove_multiple_hits ───────
OLD2 = "df[['AS']] = df[['AS']].astype('int16')"
NEW2 = (
    "df['AS'] = df['AS'].apply(\n"
    "        lambda x: int(str(x).split(':')[-1]) if ':' in str(x) else int(x)\n"
    "    ).astype('int16')"
)

if OLD2 in src:
    src = src.replace(OLD2, NEW2, 1)
    report.append("Patch 2 applied: AS tag integer extraction fixed.")
elif NEW2 in src:
    report.append("Patch 2 already applied.")
else:
    report.append("Patch 2 not needed or already handled.")

# ── Patch 3: guard in parse() before temp[5] access ───────────────────
OLD3 = '    cigr = temp[5]'
NEW3 = (
    '    if len(temp) < 11:\n'
    '        return None,None,None,None,None,None,None,None,None,0,False,False,False\n'
    '    cigr = temp[5]'
)

if OLD3 in src and NEW3 not in src:
    src = src.replace(OLD3, NEW3, 1)
    report.append("Patch 3 applied: parse() length guard added.")
elif NEW3 in src:
    report.append("Patch 3 already applied.")
else:
    report.append("Patch 3 not needed or already handled.")

# ── Patch 4: skip None returns in worker() ────────────────────────────
OLD4 = (
    '        k, r, nm, nid, rd, rs, rq, flag, cigr, read_len, '
    'pri_aln_flag, valid_match_flag, valid_acc_flag'
    ' = parse(line, matchFactor, excluded_acc_list)'
)
NEW4 = (
    '        _parsed = parse(line, matchFactor, excluded_acc_list)\n'
    '        if _parsed[0] is None:\n'
    '            continue\n'
    '        k, r, nm, nid, rd, rs, rq, flag, cigr, read_len, '
    'pri_aln_flag, valid_match_flag, valid_acc_flag = _parsed'
)

if OLD4 in src:
    src = src.replace(OLD4, NEW4, 1)
    report.append("Patch 4 applied: worker() None guard added.")
elif '_parsed = parse' in src:
    report.append("Patch 4 already applied.")
else:
    report.append("Patch 4 not needed or already handled.")

with open(path, 'w') as f:
    f.write(src)

for line in report:
    print(line)

if failed:
    for line in failed:
        print(line, file=sys.stderr)
    sys.exit(1)
PYEOF

if [ $? -ne 0 ]; then
    echo "ERROR: gottcha2.py patching failed. Printing relevant source lines for diagnosis:" >&2
    grep -n "choices\|presetx\|map-pb\|map-ont\|map-hifi" "$GOTTCHA2_BIN" | head -40 >&2 || true
    exit 1
fi

echo ""
echo "3. Running GOTTCHA2 taxonomic profiling (PacBio HiFi, direct mode)..."
echo "Command:"
echo "\"$GOTTCHA2_BIN\" \\"
echo "  -i \"$INPUT\" \\"
echo "  -d \"$DB_PREFIX\" \\"
echo "  -ti \"$DB_TAX\" \\"
echo "  -t 32 \\"
echo "  -xm map-hifi \\"
echo "  -ss 0.85,0.95 \\"
echo "  -o \"$OUTPUT_PREFIX\""
echo ""

# NOTE (2026-09): the installed GOTTCHA2 CLI has been upgraded since this
# script was written for VERSION 2.1.10. In 2.1.10, "-a <FLOAT>[,<FLOAT>]"
# set the linear-score thresholds and "-A [FILE]" was the accession
# exclusion list. The new version swapped the letters: "-a [FILE]" is now
# the accession exclusion list and "-ss <FLOAT>[,...]" carries the score
# thresholds. Using "-a 0.85,0.95" against the new CLI makes it look for a
# file literally named "0.85,0.95" and fail with
# "Accession exclusion list 0.85,0.95 not found." Fixed to "-ss 0.85,0.95".
"$GOTTCHA2_BIN" \
  -i "$INPUT" \
  -d "$DB_PREFIX" \
  -ti "$DB_TAX" \
  -t 32 \
  -xm map-hifi \
  -ss 0.85,0.95 \
  -o "$OUTPUT_PREFIX"

echo ""
echo "=========================================="
echo "GOTTCHA2 PacBio HiFi analysis completed!"
echo "End time: $(date)"
echo "=========================================="
echo ""
echo "Output files are in: ${OUTPUT_PREFIX}/"

INPUT_BASENAME=$(basename "$INPUT" .fastq.gz)
echo "Expected output files may include:"
echo "  - ${INPUT_BASENAME}.tsv        : Classification results"
echo "  - ${INPUT_BASENAME}.full.tsv   : Full results"
echo "  - ${INPUT_BASENAME}.lineage.tsv: Lineage results"
echo ""
echo "Listing output directory:"
ls -lh "${OUTPUT_PREFIX}/" || true
