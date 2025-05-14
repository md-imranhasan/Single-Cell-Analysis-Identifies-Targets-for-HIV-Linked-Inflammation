
# load into your session
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(clustree)
library(SingleR)
library(celldex)

data_dir <- 'D:/Bioinformatics/Single Cell/HIV/GSE157829_RAW/GSM4775588' # set working directory 
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx //load dataset
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix) # create seurat object
pbmc<-seurat_object
# pre-processing (start) 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells in the control group
head(pbmc@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
pbmc <- subset(pbmc, subset = nFeature_RNA >= 200)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2
#pre-processing end
# normalization
pbmc <- NormalizeData(pbmc)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
all.genes <- rownames(pbmc)
# scale data
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, verbose = FALSE)
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)
print(pbmc[["pca"]], dims = 40:50, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca",balanced=TRUE)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# removing unwanted features
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:10) 
JackStrawPlot(pbmc, dims = 1:10)

# clustering
ElbowPlot(pbmc,ndims=20) #we have to see the diagram and then change the value.. (ndims = 50)

pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE) #50 general
pbmc <- FindNeighbors(pbmc, dims = 1:20) # 50 general uses
pbmc <- FindClusters(object = pbmc,  resolution = c(0.5, 1, 1.5),  dims.use = 1:20,  save.SNN = TRUE)
#pbmc <- RunTSNE(pbmc, dims = 1:50, verbose = FALSE) # if using TSNE
#saveRDS(pbmc, file = "cluster.rds")
clustree(pbmc)
Idents(pbmc) <- pbmc$RNA_snn_res.1 # by default .1 dea ase, I can change .5 to 1.5; any poitive number
DimPlot(pbmc, reduction = "umap", label=TRUE)

all_markers = FindAllMarkers(pbmc, logfc.threshold = 0.25)
#write.csv(all_markers,"D:/imran vai/all_cluster_tsne.csv")

# if create individual files for each cluster
# ident.1="Cluster No"
markers = FindMarkers(data,ident.1 = 1, logfc.threshold = 0.25)
write.csv(markers,"D:/imran vai/cluster1.csv")

# Visualization Heatmap
all_markers  %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(
  pbmc,
  features = top10$gene,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = c("red", "blue"),
  disp.min = -2.50,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 6,
  hjust = 1.5,
  vjust = 2.5,
  angle = 30,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.01,
  combine = TRUE
) + scale_fill_gradientn(colors = c("gray", "white", "red"))

#Reference based Annotation
sce <- GetAssayData(object = pbmc, assay = "RNA", slot = "data")
refMonaco <- MonacoImmuneData()
prediction_Monaco_main <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc), labels=refMonaco$label.main)
prediction_Monaco_fine <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc), labels=refMonaco$label.fine)

predicted_Monaco <- data.frame(cluster=sort(unique(Idents(pbmc))), Monaco_main= prediction_Monaco_main$labels, Monaco_fine=prediction_Monaco_fine$labels)
predicted_Monaco
write.csv(predicted_Monaco,"D:/imran vai/ann.csv")


