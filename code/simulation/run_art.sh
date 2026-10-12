#!/bin/bash
# Call art_illumina separately for each genome; depth is given explicitly in plan_illumina.tsv
set -uo pipefail
OUT=regen/art_out; mkdir -p "$OUT"
run_one() {
  IFS=$'\t' read -r gid path cov frag sd tgt len dom <<< "$1"
  o="$OUT/${gid}"
  [ -s "${o}1.fq" ] && return 0
  art_illumina -ss HS25 -i "$path" -p -l 150 -f "$cov" -m "$frag" -s "$sd" \
               -na -rs 42 -o "$o" > "${o}.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s "${o}1.fq" ]; then
    echo -e "FAIL\t${gid}\t${dom}\t${len}\t${cov}" >> "$OUT/failures.txt"
  fi
}
export -f run_one; export OUT
: > "$OUT/failures.txt"
cat regen/art_tasks.tsv | xargs -d'\n' -P 32 -I{} bash -c 'run_one "$@"' _ {}
echo "Generation finished. Failures: $(wc -l < $OUT/failures.txt)"
