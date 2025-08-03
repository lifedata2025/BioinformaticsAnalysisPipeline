完整的单细胞RNA-seq（scRNA-seq）分析流程，适用于10x Genomics平台的数据，主要基于R语言和Seurat包实现，适合常规研究使用（比如细胞亚群识别、差异表达分析等）：

# 🧬 单细胞RNA-seq分析完整流程（基于Seurat）

作者：染山    
日期：2025年8月3日  

### 一般主要包括以下步骤：

1. 数据读取

2. 数据质控（Quality Control）

3. 数据标准化（Normalization）

4. 高变基因识别（Variable Features）

5. 数据缩放（Scaling）

6. 主成分分析（PCA）

7. 聚类分析（Clustering）

8. UMAP降维可视化

9. 差异表达分析与marker基因识别

10. 细胞类型注释（可选）

11. 可视化图形展示



---

## 🗂️ 0. 所需R包安装与加载

```r
install.packages("Seurat")
install.packages("patchwork")
install.packages("dplyr")
BiocManager::install("biomaRt")

library(Seurat)
library(patchwork)
library(dplyr)
library(biomaRt)
```


## 📁 1. 数据读取

读取10X Genomics平台的数据：

```r
data_dir <- "your_path_to/filtered_feature_bc_matrix"
sc.data <- Read10X(data.dir = data_dir)
seurat.obj <- CreateSeuratObject(counts = sc.data, project = "scRNA", min.cells = 3, min.features = 200)
```

## 🧼 2. 数据质控（Quality Control）

计算线粒体基因百分比，并进行初步筛选：

```r
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
```

## 🔬 3. 数据标准化（Normalization）

```r
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
```

## 🔍 4. 高变基因识别（Variable Features）

```r
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(seurat.obj)
```

## 📊 5. 数据缩放（Scaling）

```r
seurat.obj <- ScaleData(seurat.obj, features = rownames(seurat.obj))
```

## 📈 6. 主成分分析（PCA）

```r
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
ElbowPlot(seurat.obj)
```

## 🧱 7. 聚类分析（Clustering）

```r
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:15)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
```

## 🌐 8. UMAP降维可视化

```r
seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)
DimPlot(seurat.obj, reduction = "umap", label = TRUE)
```

## 🔎 9. 差异表达分析与marker基因识别

```r
markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
```

## 🧬 10. 细胞类型注释（可选）

可通过绘制marker基因表达图辅助注释：

```r
FeaturePlot(seurat.obj, features = c("CD3D", "MS4A1", "LYZ", "PPBP"))
```

## 🎨 11. 可视化图形展示

```r
VlnPlot(seurat.obj, features = c("CD3D"), group.by = "seurat_clusters")
DoHeatmap(seurat.obj, features = top10$gene)
```

## 💾 12. 保存分析结果

```r
saveRDS(seurat.obj, file = "seurat_obj.rds")
write.csv(markers, "cluster_markers.csv")
```
