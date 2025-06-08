#!/bin/bash

set -euo pipefail

# 配置参数
THREADS=8

if [ "$#" -lt 1 ]; then
  echo "用法: $0 SRR_ID1 [SRR_ID2 ... SRR_IDn]"
  exit 1
fi

SRA_IDS=("$@")  

mkdir -p {raw_data,ref/{hisat2_index},trimmed,alignments,counts,qc_reports/{fastqc_raw,fastqc_trimmed},results,scripts,logs}

LOG_FILE="logs/rasflow_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=============================================="
echo "Starting RASFLOW RNA-Seq Analysis Pipeline"
echo "Timestamp: $(date)"
echo "Input SRA IDs: ${SRA_IDS[@]}"
echo "=============================================="

if [ -z "${CONDITIONS+x}" ]; then
  echo " CONDITIONS变量未设置!"
  echo "在使用脚本前设置CONDITIONS数组变量，例如:"
  echo 'CONDITIONS=("control" "control" "treated" "treated")'
  echo "确保顺序与输入的SRA ID一致"
  exit 1
fi

# 步骤00: 数据下载（仅下载未处理过的样本）
echo "Step 00: Downloading SRA data..."
{
    echo "Downloading SRA data for ${SRA_IDS[@]}..."
    for SRA_ID in "${SRA_IDS[@]}"; do
      # 检查是否已处理过该样本
      if [ -f "raw_data/${SRA_ID}_1.fastq.gz" ] || [ -f "raw_data/${SRA_ID}.fastq.gz" ]; then
        echo "Sample $SRA_ID already exists. Skipping download."
        continue
      fi
      
      echo "Processing $SRA_ID ..."
      prefetch "$SRA_ID"
      fasterq-dump "$SRA_ID" -O raw_data/ --split-files --threads "$THREADS"
      gzip raw_data/"$SRA_ID"*.fastq
    done
} || {
    echo "Error in Step 00: Data download failed"
    exit 1
}

# 步骤01: 参考基因组准备（跳过已存在的）
echo "Step 01: Preparing reference genome..."
{
    cd ref
    
    # 检查参考基因组文件和索引是否存在
    if [[ -f "genome.fa" && -f "annotation.gtf" && -f "hisat2_index/genome.1.ht2" ]]; then
        echo "Reference genome and index already exist. Skipping download and index building."
        cd ..
        return 0
    fi
    
    echo "Downloading mouse reference genome..."
    wget -c ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa
    
    echo "Downloading annotation file..."
    wget -c ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
    gunzip Mus_musculus.GRCm39.110.gtf.gz
    mv Mus_musculus.GRCm39.110.gtf annotation.gtf
    
    echo "Building HISAT2 index..."
    mkdir -p hisat2_index
    hisat2-build genome.fa hisat2_index/genome
    
    cd ..
} || {
    echo "Error in Step 01: Reference preparation failed"
    exit 1
}

# 步骤02: 原始数据质控（仅对新样本）
echo "Step 02: Running FastQC on raw data..."
{
    mkdir -p qc_reports/fastqc_raw
    
    for SRA_ID in "${SRA_IDS[@]}"; do
      # 查找可能的fastq文件（单端）
      for file in raw_data/"${SRA_ID}"*.fastq.gz; do
        if [ ! -f "$file" ]; then
          continue
        fi
        
        # 检查是否已生成QC报告
        qc_report="qc_reports/fastqc_raw/$(basename "${file%.fastq.gz}")_fastqc.html"
        if [ ! -f "$qc_report" ]; then
          echo "Running FastQC on $file"
          fastqc -o qc_reports/fastqc_raw "$file"
        else
          echo "FastQC report already exists for $file. Skipping."
        fi
      done
    done
} || {
    echo "Error in Step 02: Raw data QC failed"
    exit 1
}

# 步骤03: 数据修剪（仅对新样本）
echo "Step 03: Trimming adapters with Trim Galore..."
{
    mkdir -p trimmed
    
    for SRA_ID in "${SRA_IDS[@]}"; do
      # 查找可能的fastq文件（单端）
      for file in raw_data/"${SRA_ID}"*.fastq.gz; do
        if [ ! -f "$file" ]; then
          continue
        fi
        
        # 检查是否已修剪
        trimmed_file="trimmed/$(basename "${file%.fastq.gz}")_trimmed.fq.gz"
        if [ ! -f "$trimmed_file" ]; then
          echo "Trimming $file"
          trim_galore "$file" -o trimmed --cores "$THREADS"
        else
          echo "Trimmed file already exists for $file. Skipping."
        fi
      done
    done
} || {
    echo "Error in Step 03: Trimming failed"
    exit 1
}

# 步骤04: 修剪后数据质控（仅对新样本）
echo "Step 04: Running FastQC on trimmed data..."
{
    mkdir -p qc_reports/fastqc_trimmed
    
    for file in trimmed/*_trimmed.fq.gz; do
      # 检查是否已生成QC报告
      qc_report="qc_reports/fastqc_trimmed/$(basename "${file%.fq.gz}")_fastqc.html"
      if [ ! -f "$qc_report" ]; then
        echo "Running FastQC on $file"
        fastqc -o qc_reports/fastqc_trimmed "$file"
      else
        echo "FastQC report already exists for $file. Skipping."
      fi
    done
} || {
    echo "Error in Step 04: Trimmed data QC failed"
    exit 1
}

# 步骤05: 序列比对（仅对新样本）
echo "Step 05: Aligning reads with HISAT2..."
{
    mkdir -p alignments
    
    INDEX=ref/hisat2_index/genome
    
    for SRA_ID in "${SRA_IDS[@]}"; do
      # 查找修剪后的文件
      for file in trimmed/"${SRA_ID}"*_trimmed.fq.gz; do
        if [ ! -f "$file" ]; then
          continue
        fi
        
        base=$(basename "$file" _trimmed.fq.gz)
        sam_file="alignments/${base}.sam"
        
        if [ ! -f "$sam_file" ] && [ ! -f "alignments/${base}.sorted.bam" ]; then
          echo "Aligning $file"
          hisat2 -x "$INDEX" \
              -U "$file" \
              -S "$sam_file" \
              --threads "$THREADS"
        else
          echo "Alignment already exists for $file. Skipping."
        fi
      done
    done
} || {
    echo "Error in Step 05: Alignment failed"
    exit 1
}

# 步骤06: BAM文件处理（仅对新样本）
echo "Step 06: Sorting and indexing BAM files..."
{
    for SRA_ID in "${SRA_IDS[@]}"; do
      # 查找SAM文件
      for sam in alignments/"${SRA_ID}"*.sam; do
        if [ ! -f "$sam" ]; then
          continue
        fi
        
        base=$(basename "$sam" .sam)
        bam_file="alignments/${base}.sorted.bam"
        
        if [ ! -f "$bam_file" ]; then
          echo "Processing $sam"
          samtools view -bS "$sam" | samtools sort -o "$bam_file" -@ "$THREADS"
          samtools index "$bam_file" -@ "$THREADS"
          rm "$sam"
        else
          echo "BAM file already exists for $sam. Skipping."
        fi
      done
    done
} || {
    echo "Error in Step 06: BAM processing failed"
    exit 1
}

# 步骤07: 基因计数（总是重新运行以确保包含所有样本）
echo "Step 07: Counting reads with featureCounts..."
{
    mkdir -p counts
    
    echo "Running featureCounts on all BAM files (including new ones)..."
    featureCounts -T "$THREADS" -t exon -g gene_id \
        -a ref/annotation.gtf -o counts/gene_counts.txt alignments/*.sorted.bam
    
    # 备份旧计数文件
    if [ -f "counts/gene_counts.txt" ]; then
      cp counts/gene_counts.txt "counts/gene_counts_$(date +%Y%m%d_%H%M%S).bak"
    fi
} || {
    echo "Error in Step 07: Read counting failed"
    exit 1
}

# 步骤08: 差异表达分析
echo "Step 08: Running DESeq2 for differential expression analysis..."
{
    mkdir -p results
    
    cat << 'EOF' > scripts/06_deseq2.R
library("DESeq2")
library("tximport")
library("readr")
library("ggplot2")

countData <- read.delim("counts/gene_counts.txt", comment.char="#", row.names=1)
countData <- countData[ ,6:ncol(countData)]  

conditions <- strsplit(Sys.getenv("CONDITIONS"), " ")[[1]]

colData <- data.frame(row.names=colnames(countData), 
                     condition=factor(conditions, levels=unique(conditions)))

dds <- DESeqDataSetFromMatrix(countData=countData, 
                             colData=colData, 
                             design=~condition)
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), file="results/deseq2_results.csv")

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    BiocManager::install("EnhancedVolcano")
}
library("EnhancedVolcano")
pdf("results/volcano_plot.pdf")
EnhancedVolcano(res, 
               lab=rownames(res), 
               x='log2FoldChange', 
               y='pvalue', 
               pCutoff=0.05, 
               FCcutoff=1.5,
               title='Differential Expression Analysis',
               subtitle='Volcano plot of gene expression changes')
dev.off()
EOF

    export CONDITIONS="${CONDITIONS[*]}"
    
    echo "Running DESeq2 analysis with conditions: ${CONDITIONS[@]}"
    Rscript scripts/06_deseq2.R
    
    echo "Conditions used: ${CONDITIONS[@]}" > results/conditions_used.txt
} || {
    echo "Error in Step 08: Differential expression analysis failed"
    exit 1
}

echo "=============================================="
echo "Pipeline completed successfully!"
echo "Timestamp: $(date)"
echo "Processed SRA IDs: ${SRA_IDS[@]}"
echo "Conditions used: ${CONDITIONS[@]}"
echo "=============================================="
