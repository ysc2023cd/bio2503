#!/bin/bash
set -e

INPUT_DIR=../raw_data
OUTPUT_DIR=../qc_reports/fastqc_raw
FILE_NAME="*.fastq.gz"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/$FILE_NAME; do
  fastqc -o "$OUTPUT_DIR" "$file"
done
