#!/bin/bash
set -e

INPUT_DIR=../alignments

for sam in "$INPUT_DIR"/*.sam; do
  base=$(basename "$sam" .sam)
  samtools view -bS "$sam" | samtools sort -o "$INPUT_DIR/${base}.sorted.bam"
  samtools index "$INPUT_DIR/${base}.sorted.bam"
  rm "$sam"
done

