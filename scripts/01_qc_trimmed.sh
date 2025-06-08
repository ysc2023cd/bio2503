#!/bin/bash
set -e

INPUT_DIR=../trimmed
OUTPUT_DIR=../qc_reports/fastqc_trimmed
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.fq.gz; do
  fastqc -o "$OUTPUT_DIR" "$file"
done

