# DESeq2 差异表达分析脚本
library("DESeq2")
library("tximport")
library("readr")
library("ggplot2")

countData <- read.delim("../counts/gene_counts.txt", comment.char="#", row.names=1)
countData <- countData[ ,6:ncol(countData)]  # 去掉注释列
colData <- data.frame(row.names=colnames(countData), condition=c("control", "control", "treated", "treated"))

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file="../results/deseq2_results.csv")

# 火山图示例
library("EnhancedVolcano")
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue', pCutoff=0.05, FCcutoff=1.5)
