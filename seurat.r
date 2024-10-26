
library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("results")
list.files()
list.files(path = "data/one")
getwd()
setwd("../")
sample_dir <- list.files(path = "data", pattern = "CID", full.names = TRUE)
sample_dir
aa <- tail(sample_dir, 1)
aa
cid <- list()
#Read the samples into the seurat object in a loop function
for (dir in sample_dir){
    #Read te count matrix from the directory
    #print(dir)
    project <- basename(dir)
    #print(project)
    cid_matrix = Read10X(data.dir = dir, gene.column = 1)
    cid_obj = CreateSeuratObject(counts = cid_matrix, project = project, min.cells = 3, min.features = 200)
    cid[[project]] <- cid_obj
}
cid

count_matrix_cid1 = Read10X(data.dir = "data/one/CID3586", gene.column = 1)
count_matrix_cid2 = Read10X(data.dir = "data/one/CID3838", gene.column = 1)
count_matrix_cid1
cid1 = CreateSeuratObject(counts = count_matrix_cid1, project = "CID3586", min.cells = 3, min.features = 200)
cid2 = CreateSeuratObject(counts = count_matrix_cid2, project = "CID3838", min.cells = 3, min.features = 200)

merged_cid <- merge(cid1, y = c(cid2), add.cell.ids = ls()[1:2], project = "merged_cid")
merged_cid

cell_ids <- names(cid)
merged_cid <- cid[[1]]
for (i in 2:length(cid)){
    merged_cid <- merge(merged_cid, y = cid[[i]], add.cells.ids = c(cell_ids[1], cell_ids[i], project = "merged_cid"))
}
merged_cid

merged_cid[["percent.mt"]] <- PercentageFeatureSet(merged_cid, pattern = "^MT-")
vln_plot <- VlnPlot(merged_cid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "results/mt_violin_plot.png", plot = vln_plot, width = 20, height = 6, dpi = 300)

plot1 <- FeatureScatter(merged_cid, feature1 =  "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_cid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatter_plot <- plot1 + plot2
ggsave(filename = "results/feature_scatter_plot.png", plot = scatter_plot, width = 20, height = 6, dpi = 300)

merged_cid <- subset(merged_cid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt)
merged_cid <- NormalizeData(merged_cid, normalization.method = "LogNormalize", scale.factor = 10000)
merged_cid <- NormalizeData(merged_cid)
merged_cid <- FindVariableFeatures(merged_cid, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merged_cid), 10)
plot1 <- VariableFeaturePlot(merged_cid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variable_feature_plot <- plot1 + plot2
ggsave(filename = "results/variable_feature_plot.png", plot = variable_feature_plot, width = 20, height = 6, dpi = 300)

all.genes <- rownames(merged_cid)
merged_cid <- ScaleData(merged_cid, features = all.genes)
merged_cid <- RunPCA(merged_cid, features = VariableFeatures(object = merged_cid))
print(merged_cid[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged_cid, dims = 1:2, reduction = "pca")
DimPlot(merged_cid, reduction = "pca") + NoLegend()
DimHeatmap(merged_cid, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(merged_cid)
merged_cid <- FindNeighbors(merged_cid, dims = 1:10)
merged_cid <- FindClusters(merged_cid, resolution = 0.5)
head(Idents(merged_cid), 5)
merged_cid <- RunUMAP(merged_cid, dims = 1:10)
DimPlot(merged_cid, reduction = "umap")
merged_cid <- JoinLayers(merged_cid)
cluster2.markers <- FindMarkers(merged_cid, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(merged_cid, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
merged_cid.markers <- FindAllMarkers(merged_cid, only.pos = TRUE)
merged_cid.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(merged_cid, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(merged_cid, features = c("MS4A1", "CD79A"))
