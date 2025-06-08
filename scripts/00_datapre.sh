#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
  echo "用法: $0 SRR_ID1 [SRR_ID2 ... SRR_IDn]"
  exit 1
fi

mkdir -p raw_data/

for SRA_ID in "$@"; do
  echo "处理 $SRA_ID ..."
  prefetch "$SRA_ID"
  fasterq-dump "$SRA_ID" -O raw_data/ --split-files --threads 8
  gzip raw_data/"$SRA_ID"_*.fastq
done

echo "下载完成。"


