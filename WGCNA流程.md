**WGCNA（加权基因共表达网络分析，Weighted Gene Co-expression Network Analysis）**的一个完整流程，适用于 R 环境下对转录组数据（如 RNA-seq）或微阵列数据进行共表达模块识别和关联分析。流程涵盖数据预处理、网络构建、模块识别、表型关联及功能富集分析。
🔬 WGCNA 分析完整流程（以 R 为例）
🧩 Step 1：安装并加载必要的 R 包

# 安装 WGCNA 包（只需安装一次）
install.packages("BiocManager")
BiocManager::install("WGCNA")

# 加载包
library(WGCNA)

# 为防止过多线程导致崩溃，建议关闭线程检查
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

🧬 Step 2：准备和清洗表达数据

    表达矩阵 datExpr 的格式：行为基因，列为样本。

    剔除缺失值和低表达基因。

# 载入表达数据
exprData <- read.table("expression_matrix.txt", header=TRUE, row.names=1)

# 转置矩阵，使行为样本，列为基因
datExpr0 <- as.data.frame(t(exprData))

# 检查是否有异常样本或基因
gsg <- goodSamplesGenes(datExpr0, verbose=3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

🧪 Step 3：样本聚类，检测离群样本

sampleTree <- hclust(dist(datExpr0), method="average")
plot(sampleTree, main="Sample clustering to detect outliers")

如有离群样本，可以剔除后重新赋值 datExpr
⚙️ Step 4：选择软阈值（soft-thresholding power）

powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose=5)

# 绘图查看合适 power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
     type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.9, col="blue") # 选择接近0.9的拐点

🧱 Step 5：构建网络并识别模块

softPower <- 6  # 假设你选择的是 power=6
adjacency <- adjacency(datExpr0, power=softPower)

# 将邻接矩阵转化为拓扑重叠矩阵（TOM）
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 聚类并识别模块
geneTree <- hclust(as.dist(dissTOM), method="average")
plot(geneTree, main="Gene clustering on TOM")

# 模块识别（动态剪枝）
dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM,
                             deepSplit=2, pamRespectsDendro=FALSE,
                             minClusterSize=30)
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels=FALSE)

🧠 Step 6：合并相似模块

MEList <- moduleEigengenes(datExpr0, colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method="average")

plot(METree, main="Clustering of module eigengenes")
abline(h=0.25, col="red")  # 可调节阈值

merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight=0.25, verbose=3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 更新模块颜色
moduleColors <- mergedColors

🧮 Step 7：模块与表型的相关性分析

# 读取样本表型数据
traitData <- read.table("trait_data.txt", header=TRUE, row.names=1)

# 与表达数据样本一致
datTraits <- traitData[match(rownames(datExpr0), rownames(traitData)), ]

# 计算模块与性状的皮尔逊相关
moduleTraitCor <- cor(mergedMEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))

# 可视化热图
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(mergedMEs),
               ySymbols=names(mergedMEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=signif(moduleTraitCor, 2),
               main="Module-trait relationships")

🎯 Step 8：识别核心基因（hub genes）

# 举例：假设你对turquoise模块感兴趣
module <- "turquoise"
moduleGenes <- moduleColors == module

# 计算模块成员度与性状相关性
geneModuleMembership <- cor(datExpr0, mergedMEs, use="p")
geneTraitSignificance <- cor(datExpr0, datTraits$Trait1, use="p")

# 可视化模块基因的相关性
plot(abs(geneModuleMembership[moduleGenes, module]),
     abs(geneTraitSignificance[moduleGenes]),
     xlab="Module Membership",
     ylab="Gene Trait Significance",
     main=paste("Module membership vs. trait significance\n", module))

🧬 Step 9：导出模块基因用于后续富集分析

# 导出某个模块的基因列表
turquoiseGenes <- names(datExpr0)[moduleColors == "turquoise"]
write.table(turquoiseGenes, "turquoise_genes.txt", quote=FALSE, row.names=FALSE)

之后你可以使用 clusterProfiler、DAVID、Metascape 等工具进行 GO / KEGG 富集分析。
✅ 总结：WGCNA 分析九大步骤
步骤	内容
1	安装并加载 WGCNA
2	表达数据预处理
3	样本聚类与异常检测
4	选择软阈值（scale-free）
5	构建网络，识别模块
6	合并模块
7	模块与表型相关性分析
8	寻找核心基因
9	导出模块基因并做富集分析
