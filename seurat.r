
library(Seurat)
library(dplyr)


list.files()
list.files(path = "data/one")
getwd()
setwd("../")
sample_dir <- list.files(path = "data/one", pattern = "CID", full.names = TRUE)
sample_dir
pbmc1 <- list()
#Read the samples into the seurat object in a loop function
for (dir in sample_dir){
    #Read te count matrix from the directory
    print(dir)
    count_matrix = Read10X(data.dir = dir, gene.column = 1)
    #create seurat object
    #seurat_obj <- CreateSeuratObject(counts = count_matrix, project = dir)
    seurat_obj = CreateSeuratObject(counts = count_matrix, project = dir, min.cells = 3, min.features = 200)
    #append the object to the empty list
    pbmc1[[dir]] <- seurat_obj

}

count_matrix_cid1 = Read10X(data.dir = "data/one/CID3586", gene.column = 1)
count_matrix_cid2 = Read10X(data.dir = "data/one/CID3838", gene.column = 1)

cid1 = CreateSeuratObject(counts = count_matrix_cid1, project = "CID3586", min.cells = 3, min.features = 200)
cid2 = CreateSeuratObject(counts = count_matrix_cid2, project = "CID3838", min.cells = 3, min.features = 200)

merged_cid <- merge(cid1, y = c(cid2), add.cell.ids = ls()[1:2], project = "merged_cid")
merged_cid





merged_cid[["percent.mt"]] <- PercentageFeatureSet(merged_cid, pattern = "^MT-")
VlnPlot(merged_cid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(pbmc, feature1 =  "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
