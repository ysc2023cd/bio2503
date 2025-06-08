#!/bin/bash
set -e

INPUT_DIR=../raw_data
OUTPUT_DIR=../trimmed
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.fastq.gz; do
    trim_galore "$file" -o "$OUTPUT_DIR"
done

