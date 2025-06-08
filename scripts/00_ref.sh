#!/bin/bash
# 准备小鼠参考基因组并构建HISAT2索引

cd ../ref

wget -c ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa

wget -c ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz
mv Mus_musculus.GRCm39.110.gtf annotation.gtf

mkdir -p hisat2_index
hisat2-build genome.fa hisat2_index/genome
