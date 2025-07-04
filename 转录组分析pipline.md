# 一、样本准备与测序

    RNA提取

        工具：TRIzol、RNeasy（QIAGEN）

        质量控制：Nanodrop、Qubit、Agilent Bioanalyzer（RIN ≥ 7）

    文库构建与测序

        平台：Illumina NovaSeq, HiSeq（常见格式为 paired-end，150 bp）

        输出：FASTQ 文件

# 二、RNA-Seq 数据分析流程

### 第1步：质量质控（Quality Control, QC）

- 1. 检查原始测序数据质量

    工具：

        FastQC：查看每个样本的测序质量

        MultiQC：整合多个样本的FastQC报告

- 2. 去除低质量读段与接头污染

    工具：

        Trim Galore（封装 Cutadapt + FastQC）

        fastp（速度快、自动化程度高，推荐）

输出：清洗后的 clean reads

### 第2步：序列比对（Read Mapping/Alignment）

    工具：

        HISAT2（速度快，占内存少，推荐）

        STAR（对长读段更好，但内存需求高）

输入：clean reads + reference genome

输出：比对后的 SAM / BAM 文件

### 第3步：表达定量（Expression Quantification）

    方法 1：基于比对

        工具：

            StringTie：组装转录本并计算 FPKM

            featureCounts（来自 Subread）：统计已知基因的 read count

    方法 2：无比对的伪比对方法

        工具：

            Salmon 或 Kallisto：速度快、内存低，直接输出TPM/Counts

### 第4步：差异表达分析（Differential Expression Analysis, DEA）

    输入：read count 矩阵

    工具：

        DESeq2（R包，经典选择）

        edgeR（适合低表达分析）

        limma-voom（适合重复数较少的情况）

输出：差异表达基因列表（log2FoldChange，p-value，adjusted p-value）

### 第5步：功能注释与功能富集分析（Functional Enrichment Analysis）

    GO / KEGG 富集分析

        工具（R包）：

            clusterProfiler（功能强大，支持GO/KEGG/Reactome）

            topGO（专注于GO分析）

            GSEA（Gene Set Enrichment Analysis）

    注释数据库：

        org.Hs.eg.db（人类）、org.Mm.eg.db（小鼠）等

        biomaRt（用于基因ID转换）

### 第6步：可视化分析

    常用图形：

        火山图（volcano plot）

        热图（heatmap）

        PCA图（主成分分析）

        聚类图（hierarchical clustering）

    工具：

        ggplot2（R可视化核心包）

        pheatmap 或 ComplexHeatmap（绘制热图）

        EnhancedVolcano（差异表达可视化）

        plotly（交互式可视化）

### 可选步骤：新转录本发现 / 可变剪接 / 融合基因检测

    新转录本预测：StringTie, Cufflinks

    可变剪接分析：rMATS, SUPPA2, DEXSeq

    融合基因检测：STAR-Fusion, FusionCatcher

### 常用流程化工具/框架（建议使用）

[nf-core/rnaseq](https://github.com/nf-core/RNAseq) （基于 Nextflow，社区维护的标准RNA-seq流程）

SnakeMake 或 Nextflow：可自定义的流程化执行平台

Galaxy：图形化界面操作，无需编程

**示例流程总结图（概念）**

```
FASTQ → QC (FastQC) → Trimming (fastp) → Mapping (HISAT2) →
Count (featureCounts) → DEG (DESeq2) → GO/KEGG (clusterProfiler) →
Visualization (ggplot2 / pheatmap / PCA)
```

### 验证与后续实验

    RT-qPCR验证： 选择部分关键的差异表达基因，通过**定量逆转录PCR（RT-qPCR）**进行独立的实验验证，以确认测序结果的准确性。
    其他实验验证： 根据分析结果，可以设计更深入的分子生物学实验，如Western Blot验证蛋白表达、免疫组化（IHC）或免疫荧光（IF）验证蛋白定位、基因功能研究（基因敲除/过表达）等，以进一步阐明基因的功能和机制。

