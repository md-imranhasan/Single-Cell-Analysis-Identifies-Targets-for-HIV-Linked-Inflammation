# install packages
#install.packages("Seurat")
#install.packages("Matrix")
#remove.packages("Matrix")
#install.packages("hdf5r")
#install.packages("clustree")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager",repos = "http://cran.us.r-project.org")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
# load into your session
library(dplyr)
library(Matrix)
library(Seurat)
library(patchwork)
library(clustree)
library(SingleR)
library(celldex) 
library(tidyverse)
# Read data for each donor/ sample; i have three sample

donor1<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775588")
donor2<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775589")
donor3<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775590")
donor4<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775591")
donor5<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775592")
donor6<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775593")
donor7<- Read10X(data.dir = "D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775594")
# Create Seurat objects for each  sample
seurat_donor1 <- CreateSeuratObject(counts = donor1)
seurat_donor1
seurat_donor2 <- CreateSeuratObject(counts = donor2)
seurat_donor2
seurat_donor3 <- CreateSeuratObject(counts = donor3)
seurat_donor3
seurat_donor4 <- CreateSeuratObject(counts = donor4)
seurat_donor4
seurat_donor5 <- CreateSeuratObject(counts = donor5)
seurat_donor5
seurat_donor6 <- CreateSeuratObject(counts = donor6)
seurat_donor6
seurat_donor7 <- CreateSeuratObject(counts = donor7)
seurat_donor7
# Merge the three Seurat objects
seurat_combined <- merge(seurat_donor1, y = c(seurat_donor2, seurat_donor3, seurat_donor4, seurat_donor5, seurat_donor6, seurat_donor7), add.cell.ids = ls()[1:7],  project = "seurat_combined" ) 
#joining each object into a layers
seurat_combined <- JoinLayers(seurat_combined)
seurat_combined
view(seurat_combined@meta.data)
range(seurat_combined$nFeature_RNA)
range(seurat_combined$nCount_RNA)

seurat_combined <- PercentageFeatureSet(seurat_combined, pattern = "^MT-", col.name = "percent.mt" )
view(seurat_combined@meta.data)
range(seurat_combined$percent.mt)
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat_combined
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 70 & nFeature_RNA < 2500 & nCount_RNA < 8000 & percent.mt < 10)
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) #better view
seurat_combined
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(seurat_combined)
seurat_combined <- ScaleData(seurat_combined)
all.genes <- rownames(seurat_combined)
#seurat_combined <- ScaleData(seurat_combined, features = all.genes)
seurat_combined <- RunPCA(seurat_combined, verbose = FALSE)
print(seurat_combined[["pca"]], dims = 1:10, nfeatures = 5)
print(seurat_combined[["pca"]], dims = 40:50, nfeatures = 5)
VizDimLoadings(seurat_combined, dims = 1:2, reduction = "pca",balanced=TRUE)
DimHeatmap(seurat_combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_combined, dims = 1:15, cells = 500, balanced = TRUE)
DimPlot(seurat_combined, reduction = "pca", dims = c(1,10))
seurat_combined <- JackStraw(seurat_combined, num.replicate = 100)##runing
seurat_combined <- ScoreJackStraw(seurat_combined, dims = 1:20)
JackStrawPlot(seurat_combined, dims = 1:20)
ElbowPlot(seurat_combined) ##identify top pca###
ElbowPlot(seurat_combined,ndims=20)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:20, verbose = FALSE)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)  ##dims value should be same###
seurat_combined <- FindClusters(object = seurat_combined,  resolution = c(0.5, 1, 1.5),  dims.use = 1:20,  save.SNN = TRUE)
clustree(seurat_combined)
Idents(seurat_combined) <- seurat_combined$RNA_snn_res.0.5
DimPlot(seurat_combined, reduction = "umap", label=TRUE)


cluster0.markers <- FindMarkers(seurat_combined, ident.1 = 0, logfc.threshold = 0.25) #standard###
cluster.all.markers0.5 <- FindAllMarkers(seurat_combined,assay="RNA", logfc.threshold = 0.25)
write.csv(cluster.all.markers0.5,"D:/Bioinformatics/Single Cell/HIV/allcluster.csv") ##save location###
cluster.all.markers0.5  %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5 ##cluster increase genes minimum 10

DoHeatmap(seurat_combined, features = top5$gene) + NoLegend()
sce <- GetAssayData(object = seurat_combined, assay = "RNA", slot = "data")
refMonaco <- MonacoImmuneData() ##human primary ##monaco used for blood cells
prediction_Monaco_main <- SingleR(test=sce, ref=refMonaco, clusters=Idents(seurat_combined), labels=refMonaco$label.main)
prediction_Monaco_fine <- SingleR(test=sce, ref=refMonaco, clusters=Idents(seurat_combined), labels=refMonaco$label.fine)

predicted_Monaco <- data.frame(cluster=sort(unique(Idents(seurat_combined))), Monaco_main= prediction_Monaco_main$labels, Monaco_fine= prediction_Monaco_fine$labels)
predicted_Monaco
write.csv(predicted_Monaco, "D:/Bioinformatics/Single Cell/HIV/allcluster 2.csv") #location




