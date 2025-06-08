#!/bin/bash
set -e

INPUT_DIR=../trimmed
OUTPUT_DIR=../alignments
mkdir -p "$OUTPUT_DIR"

INDEX=../ref/hisat2_index/genome

for file in "$INPUT_DIR"/*_trimmed.fq.gz; do
  base=$(basename "$file" trimmed.fq.gz)
  hisat2 -x "$INDEX" \
    -U "$file" \
    -S "$OUTPUT_DIR/${base}.sam"
done

