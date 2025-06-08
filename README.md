# RASFLOW RNA-Seq 流程分析

## 目录
- [RASFLOW RNA-Seq 流程分析](#rasflow-rna-seq-流程分析)
  - [目录](#目录)
  - [项目概述](#项目概述)
  - [项目背景与分析思路](#项目背景与分析思路)
  - [功能特点](#功能特点)
  - [分析流程步骤概览](#分析流程步骤概览)
  - [完整流程步骤说明](#完整流程步骤说明)
  - [快速开始](#快速开始)
    - [系统要求](#系统要求)
  - [使用说明](#使用说明)
    - [必需参数](#必需参数)
    - [条件变量示例](#条件变量示例)
    - [使用示例](#使用示例)
      - [首次运行（两个样本）：](#首次运行两个样本)
      - [添加新样本（自动合并分析）：](#添加新样本自动合并分析)
  - [目录结构](#目录结构)
  - [注意事项](#注意事项)
  - [结果分析](#结果分析)
  - [反思与改进](#反思与改进)
  - [参考文献](#参考文献)
  - [贡献](#贡献)
  - [许可证](#许可证)
  - [github](#github)
## 项目概述

**RASFLOW** 是一个全自动化的 RNA-Seq 分析流程脚本，从原始 SRA 数据下载开始，到差异表达分析结束，提供完整的分析解决方案。支持增量添加新样本并重新分析。

---

## 项目背景与分析思路
经典型Ehlers-Danlos综合征(cEDS)是结缔组织与肌肉骨骼系统常见的遗传性疾病，其中一个变异的位点是V型胶原蛋白基因COL5A1。其一种严重的症状是伤口愈合缺陷。本组查询文献，找到一个项目，用小鼠构建了cEDS模型，用于研究伤口愈合缺陷的病理机制。本项目使用小组开发的RNA-Seq数据分析脚本，根据此实验上传在公共数据库上的测序结果进行RNA-Seq分析，尝试对野生型小鼠(WT)与模型小鼠(CKO)小鼠在受伤后同一时间的全部RNA表达情况进行分析比对，验证其模型成功构建。

## 功能特点

-  **全自动化流程**：从数据下载到差异分析一键完成
-  **增量处理**：智能跳过已处理样本，仅处理新数据
-  **参考基因组自动管理**：自动下载并构建索引
-  **质量控制**：提供原始数据与修剪后数据的 FastQC 报告
-  **差异表达分析**：基于 DESeq2 进行统计分析
-  **完整日志记录**：详细记录每一步处理信息，便于追踪和调试

---

## 分析流程步骤概览

-  **数据下载**：从 NCBI SRA 获取原始数据（fasterq-dump）
-  **参考基因组准备**：自动下载 GRCm39 基因组与注释并构建 HISAT2 索引
-  **质量控制**：原始与修剪后数据分别使用 FastQC 检查
-  **数据修剪**：使用 Trim Galore 去除低质量序列和接头污染
-  **比对分析**：HISAT2 将清洗后的序列比对至参考基因组
-  **表达量统计**：用 featureCounts 生成每个基因的 read count 表
-  **差异表达分析**：调用 R 脚本运行 DESeq2 分析表达变化，生成 CSV 文件和火山图

---

## 完整流程步骤说明

每部分的bash脚本封装在了scripts文件夹中，但是数据下载时的参数输入与完整的自动化脚本稍有不同。

1. **数据下载(00_datapre.sh)**
- 从NCBI SRA下载原始测序数据
- 支持多SRA ID输入
- 自动转为fastq.gz格式
```bash
if [ "$#" -lt 1 ]; then
  echo "用法: $0 SRR_ID1 [SRR_ID2 ... SRR_IDn]"
  exit 1
fi
```

2. **参考基因组准备(00_ref.sh)**
- 自动下载小鼠参考基因组(GRCm39)
- 下载基因注释文件(GTF)
- 构建HISAT2索引

3. **数据质控(01_qc.sh)**
- 使用FastQC进行原始数据质量评估
- 生成HTML格式的质量报告

4. **数据修剪(02_trim_single.sh)**
- 使用Trim Galore自动修剪:
- 去除adapter序列
- 去除低质量reads
- 自动检测质量编码

5. **修剪后质控(01_qc_trimmed.sh)**
- 对修剪后的数据再次进行FastQC分析
- 比较修剪前后的质量变化

6. **序列比对(03_align.sh)**
- 使用HISAT2将reads比对到参考基因组
- 输出SAM格式比对结果

7. **BAM文件处理(04_bam.sh)**
- 转换SAM为BAM格式
- 按坐标排序
- 建立索引
- 自动清理中间文件

8. **基因计数(05_count.sh)**
- 使用featureCounts进行基因水平定量
- 基于exon区域计数
- 输出基因计数矩阵

9. **差异表达分析(06_deseq2.sh&06_deseq2.R)**
- 使用DESeq2进行差异分析
- 自动生成:
  - 差异表达结果表(CSV)
  - 火山图(PDF)
  - MA图(PDF)

## 快速开始

### 系统要求

- Linux 环境
- 需安装并配置好以下软件（加入 PATH）：
  - SRA Toolkit (`prefetch`, `fasterq-dump`)
  - HISAT2
  - samtools
  - Trim Galore
  - FastQC
  - featureCounts
  - R（含 `DESeq2` 与 `EnhancedVolcano` 包）

- 建议使用conda创建虚拟环境:
  
```bash
git clone git@github.com:ysc2023cd/bio2503.git

conda env create -f  ras_env.yml
conda activate ras_env
```

## 使用说明

### 必需参数

- 至少提供一个 SRA ID 作为输入参数
- 必须事先设置 `CONDITIONS` 环境变量，顺序需与输入 SRA ID 完全对应

### 条件变量示例

```bash
CONDITIONS=("group1" "group1" "group2" "group2")
```

### 使用示例

#### 首次运行（两个样本）：

```bash
CONDITIONS=("control" "treated") ./rasflow_pipeline.sh SRR123456 SRR789012
```

#### 添加新样本（自动合并分析）：

```bash
CONDITIONS=("control" "treated" "treated") ./rasflow_pipeline.sh SRR123456 SRR789012 SRR345678
```

## 目录结构

```text
rasflow/
├── raw_data/           # 原始测序数据（fastq.gz）
├── ref/                # 参考基因组与注释文件、HISAT2 索引
│   └── hisat2_index/
├── trimmed/            # Trim Galore 处理后的数据
├── alignments/         # HISAT2 比对结果（BAM文件）
├── counts/             # featureCounts 生成的计数矩阵
├── qc_reports/         # FastQC 报告
│   ├── fastqc_raw/     # 原始数据 QC
│   └── fastqc_trimmed/ # 修剪后数据 QC
├── results/            # 分析输出
│   ├── deseq2_results.csv  # 差异表达分析结果
│   └── volcano_plot.pdf    # 火山图可视化
├── scripts/            # 自动生成的 R 脚本
└── logs/               # 运行日志文件
```

## 注意事项

-  条件变量必须在运行前显式设置，并与 SRA ID 顺序一致
-  已处理样本将跳过重复分析，仅对新增样本执行全流程
-  差异表达分析和计数矩阵总是会重算以保证样本完整性
-  大型基因组或多样本情况下，请预留足够内存和磁盘空间

## 结果分析
差异表达分析结果表明，COLV基因在单个COLV CKO样本中对比WT野生型的表达显著下降，说明模型小鼠构建成功。

## 反思与改进
-  仅选取了一对样本，分析得到的结果没有统计意义，只能初步说明模型构建成功，应该选取全部的10（COLV CKO）+7（WT）组样本。
-  可以进一步处理原实验所做的野生型成纤维细胞注射治疗或cilengitide处理后的表达情况与模型与野生型对比，观察各自对于伤口愈合缺陷的治疗效果。
-  对于脚本本身，可以加入对双端测序数据的分析。

## 参考文献
[1]Kelly-Scumpia KM, Archang MM, Purbey PK, Yokota T, Wu R, McCourt J, Li S, Crosbie RH, Scumpia PO, Deb A. Modulating the extracellular matrix to treat wound healing defects in Ehlers-Danlos syndrome. iScience. 2024 Aug 6;27(9):110676. doi: 10.1016/j.isci.2024.110676. PMID: 39262784; PMCID: PMC11389543.
[2]https://github.com/zhxiaokang/RASflow

## 贡献
本实验过程和报告由易思成、方淏熠共同完成。

## 许可证
本项目采用 MIT 许可证 - 详情请见 LICENSE 文件

## [github](https://github.com/ysc2023cd/bio2503.git "跳转github")
https://github.com/ysc2023cd/bio2503.git
