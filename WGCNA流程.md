WGCNA（加权基因共表达网络分析，Weighted Gene Co-expression Network Analysis）的一个完整流程，适用于 R 环境下对转录组数据（如 RNA-seq）或微阵列数据进行共表达模块识别和关联分析。
流程涵盖数据预处理、网络构建、模块识别、表型关联及功能富集等分析。

WGCNA 分析构建的是加权的无尺度网络，通过加权的方法对相关性值进行幂次运算，这种处理方式强化了强相关，弱化了弱相关或负相关，使得相关性数值更符合无标度网络特征，更具有生物意义。
具体请查找一些关于 WGCNA 的原理叙述。

# WGCNA 分析完整流程（以 R 为例）

## Step 1：安装并加载必要的 R 包

### 安装 WGCNA 包（只需安装一次）

```
install.packages("BiocManager")
BiocManager::install("WGCNA")
```

### 加载包

```
library(WGCNA)
```

### 为防止过多线程导致崩溃，建议关闭线程检查

```
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
```

## Step 2：准备和清洗表达数据

    表达矩阵 datExpr 的格式：行为基因，列为样本。

    剔除缺失值和低表达基因。

### 载入表达数据

```
exprData <- read.table("expression_matrix.txt", header=TRUE, row.names=1)
```

### 转置矩阵，使行为样本，列为基因

```
datExpr0 <- as.data.frame(t(exprData))
```

### 检查是否有异常样本或基因

```
gsg <- goodSamplesGenes(datExpr0, verbose=3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}
```

## Step 3：样本聚类，检测离群样本

```
sampleTree <- hclust(dist(datExpr0), method="average")
plot(sampleTree, main="Sample clustering to detect outliers")

如有离群样本，可以剔除后重新赋值 datExpr
```

## Step 4：选择软阈值（soft-thresholding power）

```
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose=5)

# 绘图查看合适 power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
     type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.9, col="blue") # 选择接近0.9的拐点
```

## Step 5：构建网络并识别模块

```
softPower <- 6  # 假设你选择的是 power=6
adjacency <- adjacency(datExpr0, power=softPower)
```

### 将邻接矩阵转化为拓扑重叠矩阵（TOM）

```
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
```

### 聚类并识别模块

```
geneTree <- hclust(as.dist(dissTOM), method="average")
plot(geneTree, main="Gene clustering on TOM")
```

### 模块识别（动态剪枝）

```
dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM,
                             deepSplit=2, pamRespectsDendro=FALSE,
                             minClusterSize=30)
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels=FALSE)
```

## Step 6：合并相似模块

```
MEList <- moduleEigengenes(datExpr0, colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method="average")

plot(METree, main="Clustering of module eigengenes")
abline(h=0.25, col="red")  # 可调节阈值

merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight=0.25, verbose=3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
```

### 更新模块颜色

```
moduleColors <- mergedColors
```

## Step 7：模块与表型的相关性分析

### 读取样本表型数据

```
traitData <- read.table("trait_data.txt", header=TRUE, row.names=1)
```

### 与表达数据样本一致

```
datTraits <- traitData[match(rownames(datExpr0), rownames(traitData)), ]
```

### 计算模块与性状的皮尔逊相关

```
moduleTraitCor <- cor(mergedMEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))
```

### 可视化热图

```
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(mergedMEs),
               ySymbols=names(mergedMEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=signif(moduleTraitCor, 2),
               main="Module-trait relationships")
```

## Step 8：识别核心基因（hub genes）

### 举例：假设你对turquoise模块感兴趣

```
module <- "turquoise"
moduleGenes <- moduleColors == module
```

### 计算模块成员度与性状相关性

```
geneModuleMembership <- cor(datExpr0, mergedMEs, use="p")
geneTraitSignificance <- cor(datExpr0, datTraits$Trait1, use="p")
```

### 可视化模块基因的相关性

```
plot(abs(geneModuleMembership[moduleGenes, module]),
     abs(geneTraitSignificance[moduleGenes]),
     xlab="Module Membership",
     ylab="Gene Trait Significance",
     main=paste("Module membership vs. trait significance\n", module))
```

## Step 9：导出模块基因用于后续富集分析

### 导出某个模块的基因列表

```
turquoiseGenes <- names(datExpr0)[moduleColors == "turquoise"]
write.table(turquoiseGenes, "turquoise_genes.txt", quote=FALSE, row.names=FALSE)
```

## 总结：WGCNA 分析九大步骤

| 步骤 | 内容 |
|------|------|
|1 | 安装并加载 WGCNA |
|2 | 表达数据预处理 |
|3 | 样本聚类与异常检测 |
|4 | 选择软阈值（scale-free） |
|5 | 构建网络，识别模块 |
|6 | 合并模块 |
|7 | 模块与表型相关性分析 |
|8 | 寻找核心基因 |
|9 | 导出模块基因并做富集分析 |









    expression_matrix.txt：行为基因，列为样本（未转置）。

    trait_data.txt：行为样本，列为表型（如病理分型、分组信息等）。

你可直接运行该脚本进行WGCNA分析。
🧬 WGCNA 分析全流程模板脚本（含图示）

# -----------------------------
# 0. 加载环境
# -----------------------------
options(stringsAsFactors = FALSE)
library(WGCNA)
allowWGCNAThreads()

# -----------------------------
# 1. 载入表达数据并预处理
# -----------------------------
# 表达矩阵行为基因，列为样本
exprData <- read.table("expression_matrix.txt", header=TRUE, row.names=1)
datExpr0 <- as.data.frame(t(exprData))

# 检查缺失值和低质量样本
gsg <- goodSamplesGenes(datExpr0, verbose=3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# -----------------------------
# 2. 样本聚类检测离群值
# -----------------------------
sampleTree <- hclust(dist(datExpr0), method="average")
png("01_sample_clustering.png", width=800, height=600)
plot(sampleTree, main="Sample clustering", xlab="", sub="")
abline(h=120, col="red")  # 可根据聚类树选择截断高度
dev.off()

# -----------------------------
# 3. 选择合适的软阈值
# -----------------------------
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector=powers, verbose=5)

# 绘制 scale-free 拟合度图
png("02_scale_independence.png", width=800, height=600)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
     type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.9, col="blue")
dev.off()

# -----------------------------
# 4. 网络构建与模块识别
# -----------------------------
softPower <- 6  # 手动设置选定的power
adjacency <- adjacency(datExpr0, power=softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method="average")
dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM,
                             deepSplit=2, pamRespectsDendro=FALSE,
                             minClusterSize=30)
dynamicColors <- labels2colors(dynamicMods)

# 模块可视化
png("03_gene_dendrogram_modules.png", width=800, height=600)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels=FALSE)
dev.off()

# -----------------------------
# 5. 合并相近模块
# -----------------------------
MEList <- moduleEigengenes(datExpr0, colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method="average")

png("04_module_eigengene_clustering.png", width=800, height=600)
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")
abline(h=0.25, col="red")
dev.off()

merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight=0.25, verbose=3)
moduleColors <- merge$colors
mergedMEs <- merge$newMEs

# 更新模块可视图
png("05_merged_modules_dendrogram.png", width=800, height=600)
plotDendroAndColors(geneTree, moduleColors, "Merged module colors", dendroLabels=FALSE)
dev.off()

# -----------------------------
# 6. 与表型数据关联
# -----------------------------
traitData <- read.table("trait_data.txt", header=TRUE, row.names=1)
datTraits <- traitData[match(rownames(datExpr0), rownames(traitData)), ]

moduleTraitCor <- cor(mergedMEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))

# 热图展示模块与性状的相关性
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)

png("06_module_trait_heatmap.png", width=1000, height=800)
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(mergedMEs),
               ySymbols=names(mergedMEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               main="Module-trait relationships")
dev.off()

# -----------------------------
# 7. 提取特定模块核心基因
# -----------------------------
module <- "turquoise"  # 替换为你感兴趣的模块
moduleGenes <- moduleColors == module

geneModuleMembership <- cor(datExpr0, mergedMEs, use="p")
geneTraitSignificance <- cor(datExpr0, datTraits$Trait1, use="p")  # 修改 Trait1 为实际列名

# 可视化模块核心性
png("07_module_vs_trait_scatter.png", width=800, height=600)
plot(abs(geneModuleMembership[moduleGenes, module]),
     abs(geneTraitSignificance[moduleGenes]),
     xlab="Module Membership", ylab="Gene Significance",
     main=paste("Module membership vs. trait significance\n", module),
     col=module)
dev.off()

# -----------------------------
# 8. 导出模块基因列表
# -----------------------------
moduleGeneList <- names(datExpr0)[moduleColors == module]
write.table(moduleGeneList, paste0(module, "_genes.txt"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

📁 生成图文件说明
文件名	含义
01_sample_clustering.png	样本聚类图，检查离群
02_scale_independence.png	选择软阈值图
03_gene_dendrogram_modules.png	初步模块树状图
04_module_eigengene_clustering.png	模块特征向量聚类图
05_merged_modules_dendrogram.png	合并后的模块图
06_module_trait_heatmap.png	模块与表型的相关性热图
07_module_vs_trait_scatter.png	特定模块的hub gene散点图
✅ 使用提示

    使用时请确保你的 expression_matrix.txt 和 trait_data.txt 的行名能一一匹配。

    如果你有多个性状变量（如不同的疾病状态/分组），可以做多次 geneTraitSignificance 计算。

    可扩展用于富集分析（建议使用 clusterProfiler 包）。
