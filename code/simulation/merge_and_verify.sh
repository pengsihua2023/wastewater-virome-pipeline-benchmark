#!/bin/bash
# Merge per-genome FASTQ files and run the per-member acceptance check (guide step 7, gate 1)
set -euo pipefail
OUT=regen/art_out
S=regen/sample_708R_0     # R = regenerated
echo "[1/3] merging R1/R2 ..."
cat $OUT/genome_*1.fq | gzip -1 > ${S}_01.fq.gz &
cat $OUT/genome_*2.fq | gzip -1 > ${S}_02.fq.gz &
wait
echo "[2/3] counting reads per source accession ..."
zcat ${S}_01.fq.gz \
  | awk 'NR%4==1{h=substr($1,2); n=index(h,"-"); print (n>0 ? substr(h,1,n-1) : h)}' \
  | sort | uniq -c | sort -rn > regen/acc_counts_new.txt
echo "[3/3] done"
ls -l ${S}_0*.fq.gz | awk '{printf "  %.2f GB  %s\n", $5/1073741824, $9}'
