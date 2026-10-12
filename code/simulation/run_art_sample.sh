#!/bin/bash
# Usage: run_art_sample.sh <sample number> <random seed>
set -uo pipefail
S=$1; SEED=$2
OUT=regen/art_out_s${S}; mkdir -p "$OUT"
run_one() {
  IFS=$'\t' read -r gid path cov frag sd tgt len dom <<< "$1"
  o="$OUT/${gid}"
  [ -s "${o}1.fq" ] && return 0
  art_illumina -ss HS25 -i "$path" -p -l 150 -f "$cov" -m "$frag" -s "$sd" \
               -na -rs "$SEED" -o "$o" > "${o}.log" 2>&1
  if [ $? -ne 0 ] || [ ! -s "${o}1.fq" ]; then
    echo -e "FAIL\t${gid}\t${dom}\t${len}\t${cov}" >> "$OUT/failures.txt"
  fi
}
export -f run_one; export OUT SEED
: > "$OUT/failures.txt"
cat regen/art_tasks_s${S}.tsv | xargs -d'\n' -P 32 -I{} bash -c 'run_one "$@"' _ {}
echo "sample_${S} generation finished, failures: $(wc -l < $OUT/failures.txt)"
