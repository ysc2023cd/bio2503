#!/bin/bash
set -e

INPUT_DIR=../alignments
OUTPUT_DIR=../counts
GTF=../ref/annotation.gtf
mkdir -p "$OUTPUT_DIR"

featureCounts -T 4 -t exon -g gene_id \
  -a "$GTF" -o "$OUTPUT_DIR/gene_counts.txt" "$INPUT_DIR"/*.sorted.bam

