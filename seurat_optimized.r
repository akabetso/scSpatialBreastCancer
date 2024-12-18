library(Seurat)
library(sctransform)
library(dplyr)
library(ggplot2)
library(glmGamPoi)
library(gtable)
library(multtest)
library(metap)
library(gridExtra)
library(harmony)
library(grid)
library(patchwork)
library(infercnv)
library(future)
library(hdf5r)
library(Matrix)
library(rhdf5)
library(DropletUtils)
library(ggplot2)


#process and scale data at once
dir.create("subset_3")
dir.create("data/subset_3")
sample_dir <- list.files(path = "data", pattern = "CID", full.names = TRUE)
sample_dir
cid <- list()
#Read the samples into the seurat object in a loop function
for (dir in sample_dir){
    project <- basename(dir)
    cid_matrix <- ReadMtx(
        mtx = paste0(dir, "/matrix.mtx.gz"), 
        features = paste0(dir, "/features.tsv.gz"), 
        cells = paste0(dir, "/barcodes.tsv.gz"), 
        feature.column = 1)
    cid_obj = CreateSeuratObject(counts = cid_matrix, min.cells = 3, min.features = 250)
    cid[[project]] <- cid_obj
}
merged_cid <- merge(cid[[1]], 
                    y = cid[-1],
                    project = "merged_cid")

merged_cid[["percent.mt"]] <- PercentageFeatureSet(merged_cid, pattern = "^MT-")
quality_control_one <- VlnPlot(merged_cid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "subset_3/quality_control.png", plot = quality_control_one, width = 20, height = 6, dpi = 300)
filtered_cid <- subset(merged_cid, subset = nFeature_RNA > 200 & nCount_RNA > 250 & percent.mt < 20)

#let's rather use sctransform for normalization
#filtered_cid <- NormalizeData(filtered_cid)
#filtered_cid <- FindVariableFeatures(filtered_cid, selection.method = "vst", nfeatures = 2000) 
#filtered_cid <- ScaleData(filtered_cid) 
#top10_vf <- head(VariableFeatures(filtered_cid), 10)
#plot1 <- VariableFeaturePlot(filtered_cid)
#label_plot1 <- LabelPoints(plot = plot1, points = top10_vf, repel = TRUE)
#png("results/unintegrated/top_variable_genes.png")
#print(label_plot1)
#dev.off()

#Use alternatives for nomrlaization and scalling: SCTransform
#computationally expensive.
#sct_cid <- subset(merged_cid, subset = nFeature_RNA > 200 & nCount_RNA > 250 & percent.mt < 10)
sct_cid <- SCTransform(filtered_cid, vars.to.regress = "percent.mt", variable.features.n = 2000)
#sct_top10_vf <- head(VariableFeatures(sct_cid), 10)
#plot1 <- VariableFeaturePlot(sct_cid)
#label_plot1 <- LabelPoints(plot = plot1, points = sct_top10_vf, repel = TRUE)
#png("results/unintegrated/sct_top_variable_genes.png")
#print(label_plot1)
#dev.off()

#Linear dimension reduction pca
sct_cid <- RunPCA(sct_cid, npcs = 30, verbose = FALSE)
sct_cid@reductions
pca_plot <- DimPlot(sct_cid, reduction = "pca", group.by="orig.ident")
ggsave(filename = "subset_3/unintegrated/pca.png", plot = pca_plot)
elbow_plot <- ElbowPlot(sct_cid)
png("subset_3/unintegrated/elbow_plot.png")
print(elbow_plot)
dev.off()

#dimensonal reduction owith umap
sct_cid <- RunUMAP(sct_cid, reductions = "pca", dims = 1:20)
umap_plot <- DimPlot(sct_cid, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "subset_3/unintegrated/umap_plot.png", plot = umap_plot)

#Let's find cluster or communities
#cluster_cid <- FindNeighbors(filtered_cid, dims = 1:15)
#cluster_cid <- FindClusters(cluster_cid, resolution = 0.4)
#find tha batch effect on clusters check biological to technical variance
#table(Cluster = cluster_cid$seurat_clusters, Batch = cluster_cid$orig.ident)
#batch_table <- table(Cluster = cluster_cid$seurat_clusters, Batch = cluster_cid$orig.ident)
#batch_table_top10 <- batch_table[, 1:8]
#png("results/unintegrated/batch_tabulation_table.png")
#library(grid)
#grid.table(batch_table_top10)
#dev.off()

#The tabulation presents technical variance over biological variance as in one cluster, cells in one batch significantly varies than others.
#we introduce batch effect to correct this by integration.
#loop through the samples and normalize high variables individually
#pre_integrated_cid <- lapply(c(cid[-1]), function(x) {
#  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})
#select features that are repeatedly variable accross dataset for integration
#features <- SelectIntegrationFeatures(object.list = pre_integrated_cid, nfeatures=10000)
#perform integration
#cid_anchors <- FindIntegrationAnchors(object.list = pre_integrated_cid, anchor.features = features)
#integrated_cid <- IntegrateData(anchorset = cid_anchors)

#Batch correction
#The CCA method above  is computationally expensive, so we try two more alternatives below
integrated_cid <- lapply(X = cid, FUN = function(x) {
    x <- SCTransform(x, verbose = TRUE, variable.features.n = 2000)
})
cid_features <- SelectIntegrationFeatures(object.list =  integrated_cid, nfeatures = 2000)
cid_anchors <- FindIntegrationAnchors(object.list = integrated_cid, anchor.features = cid_features)
cid_integrated <- IntegrateData(anchorset = cid_anchors, normalization.method = "SCT")
cid_integrated <- ScaleData(cid_integrated, verbose = FALSE)
cid_integrated <- RunPCA(cid_integrated, npcs = 30, verbose = FALSE)
pca_plot_intd <- DimPlot(cid_integrated, reduction = "pca", group.by = "orig.ident")
ggsave(filename = "subset_3/integrated/pca.png", plot = pca_plot_intd)
elbow_plot <- ElbowPlot(cid_integrated)
png("subset_3/integrated/elbow_plot.png")
print(elbow_plot)
dev.off()
cid_integrated <- RunUMAP(cid_integrated, reduction = "pca", dims = 1:18, n.neighbors = 10)
umap_plot <- DimPlot(cid_integrated, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "subset_3/integrated/umap_plot.png", plot = umap_plot)

#Harmony
#harmony_cid <- IntegrateLayers(object = filtered_cid, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', verbose = FALSE)
#harmony_cid@assays
#save <- harmony_cid
#harmony_cid <- RunUMAP(harmony_cid, reduction = "harmony", dims = 1:15)
#umap_harmony_plot <- DimPlot(harmony_cid, group.by = "orig.ident")
#ggsave(filename = "results/integrated/umap_plot.png", plot = umap_harmony_plot)
#perform clustering
#cluster_harmony_cid <- FindNeighbors(harmony_cid, reduction = "harmony", dims = 1:15) # Adjust dims as needed
#cluster_harmony_cid <- FindClusters(cluster_harmony_cid, resolution = 0.4) # Adjust resolution as needed
#batch_tabulation <- table(Cluster = cluster_harmony_cid$seurat_clusters, Batch = cluster_harmony_cid$orig.ident)
#batch_table_top10 <- batch_tabulation[, 1:8]
#png("results/integrated/batch_tabulation_table.png")
#library(grid)
#grid.table(batch_table_top10)
#dev.off()
#cluster_plot <- DimPlot(cluster_harmony_cid, reduction = "umap", label = TRUE)
#ggsave(filename = "results/integrated/cluster_plot.png", cluster_plot)


#Identify cell types
DefaultAssay(cid_integrated) <- "integrated"
cid_integrated <- FindNeighbors(cid_integrated, dims = 1:10)
cid_integrated <- FindClusters(cid_integrated, resolution = 0.4, algorithm = 2)
umap_plot <- DimPlot(cid_integrated, reduction = "umap", label = TRUE)
ggsave(filename = "subset_3/integrated/cluster_louvain2_plot.png", plot = umap_plot)


#plot marker genes as described in the article
DefaultAssay(cid_integrated) <- "SCT"
features <- c("EPCAM", "MKI67", "CD3D", "CD68", "MS4A1", "JCHAIN", "PECAM1", "PDGFRB")
marker_plot <- FeaturePlot(cid_integrated, features = features, pt.size = 0.1, label = TRUE)
ggsave(filename = "subset_3/integrated/marker_expression_plot.png", plot = marker_plot)
marker_vln_plot <- VlnPlot(cid_integrated, features = features)
ggsave(filename = "subset_3/integrated/marker_vln_plot.png", plot = marker_vln_plot)
cid_integrated <- PrepSCTFindMarkers(cid_integrated)
c3_markers <- FindMarkers(cid_integrated, ident.1 = 3)
head(c3_markers)
c3_marker_plot <- FeaturePlot(cid_integrated, reduction = "umap", features = rownames(head(c3_markers)), pt.size = 0.1)
ggsave(filename = "subset_3/integrated/c3_marker_expression_plot.png", plot = c3_marker_plot)
all_markers <- FindAllMarkers(cid_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC) %>% pull(gene)
top_markers_plot <- FeaturePlot(cid_integrated, features = top.markers, reduction = "umap")
ggsave(filename = "subset_3/integrated/top_marker_cluster.png", plot = top_markers_plot)
new_cluster_ids <- c("T-cells", "Dendric cells", "Endometrial stromal cells", "NK-cells", "Basal cells",
                     "Macrophages", "Fibroblast", "smooth muscles", "Endothelial cells", "B-cells", "Plasma cells", "Platelates")
names(new_cluster_ids) <- levels(cid_integrated)
top_marker_vln_plot <- VlnPlot(cid_integrated, features = top.markers)
ggsave(filename = "subset_3/integrated/top_marker_vln_plot.png", plot = top_marker_vln_plot)
#names(new_cluster_ids) <- NULL
#cid_integrated <- NULL
cid_integrated <- RenameIdents(cid_integrated, new_cluster_ids)
cluster_annotation <- DimPlot(cid_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "subset_3/integrated/cluster_annotation_plot.png", plot = cluster_annotation)
batch_tabulation <- table(Cluster = cid_integrated@active.ident, Batch = cid_integrated$orig.ident)
png("subset_3/integrated/batch_tabulation_table.png")
grid.table(batch_tabulation)
dev.off()
#DEG
cid_integrated[["cell_labels"]] <- cid_integrated@active.ident
cid_3586_vs_3946 <- FindMarkers(cid_integrated, assay = "SCT", ident.1 = "CID3586", ident.5 = "CID3946", group.by = "orig.ident", min.pct = 0.5)
head(cid_3586_vs_3946)
cid_integrated <- SetIdent(cid_integrated, value = "orig.ident")
f1 <- FeaturePlot(subset(cid_integrated, idents = c("CID3586")), features = c("CXCL13")) + ggtitle("Sample: CID3586, Marker: CXCL13")
f2 <- FeaturePlot(subset(cid_integrated, idents = c("CID3946")), features = c("CXCL13")) + ggtitle("Sample: CID3946, Marker: CXCL13")
plot <- f1 + f2
ggsave(filename = "subset_3/integrated/endometrial_CID_3586_vs_3946.png", plot = plot)


#annotate inferCNV file text
threshold <- 0.5
#infercnv_file <- data.frame(cell = colnames(cid_integrated), group = c("Unknown", "Immune"))
immune_genes <- c("MS4A1", "VPREB3", "CD3D", "GNLY", "TNF", "CCR7", "CD68")
DefaultAssay(cid_integrated) <- "RNA"
cid_integrated$immune <- apply(cid_integrated@assays$RNA$counts[immune_genes, ], 2, function(x) any(x > threshold)) #get features counts
cid_integrated$annotaion <- ifelse(cid_integrated$immune == TRUE, "Immune", "Unknown") #define cells (logical)
infer_file <- data.frame(cell = colnames(cid_integrated), group = cid_integrated$annotaion) #make df and save 
write.table(infer_file, file = "infercnv_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#get expresion matrix for infercnv
expression_matrix <- GetAssayData(cid_integrated, assay = "RNA", layer = "counts")
write.table(as.matrix(expression_matrix), file = "expresion_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) #save matrix subset
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="data/expresion_matrix.txt",
                                    annotations_file="data/infercnv_file.txt",
                                    delim="\t",
                                    gene_order_file="data/hg38_gencode_v27.txt",
                                    ref_group_names=c("Immune"))
# perform infercnv operations to reveal cnv signal
dir.create("data/subset_5")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = "data/subset_3", 
                             cluster_by_groups = T,   # cluster
                             denoise = T,
                             HMM = T
                             ) 
# plot cnv
#infercnv_plot <- plot_cnv(infercnv_obj, out_dir = "data/subset_3", obs_title = "Observation (Cells)", ref_title = "Reference (Cells)", cluster_by_groups = TRUE,
#   x.center = 1, x.range = "auto", hclust_method = 'ward.D', color_safe_pal = FALSE, output_filename = "infercnv_plot", output_format = "png", png_res = 300, dynamic_resize = 0)
infercnv_cell_plot <- DimPlot(cid_integrated, reduction = "umap", group.by = "annotaion", label = TRUE) + ggtitle("cluster by inferCNV") + theme_minimal()
ggsave(filename = "subset_3/integrated/infercnv_celltype_plot.png", plot = infercnv_cell_plot)

# Garnett Classifier
# Extraact pretrain cell type from xCell2.
library(xCell2)
library(garnett)
#library(monocle)
library(monocle3)
library(SeuratWrappers)
library(org.Hs.eg.db)
library(tidyr)
dir.create("garnett_classifier")
ref_url <- "https://dviraran.github.io/xCell2refs/references/PanCancer.xCell2Ref.rds" #download pre-trained ct
local_filename <- "garnett_classifier/PanCancer.xCell2Ref.rds"
download.file(ref_url, local_filename, mode = "wb")
PanCancer.xCell2Ref <- readRDS(local_filename)
slotNames(PanCancer.xCell2Ref)
# We have classifier from cell, let's proceed with garnett
classifier <- readRDS(local_filename)
signatures <- PanCancer.xCell2Ref@signatures
# convert seurat object to monocle celldataset
#cid_cds <- importCDS(cid_integrated, import_all = TRUE)
#cid_cds <- as.cell_data_set(cid_integrated) #convert to monocle cds
expm <- cid_integrated@assays$RNA$counts
#cells <- cid_integrated@meta.data
#genes <- rownames(expm)
#gene_metadata <- data.frame(gene_short_name = genes)
#rownames(gene_metadata) <- genes
#genes <- cid_integrated@assays$RNA$counts@Dimnames[[1]]
#cid_cds <- new_cell_data_set(expm, cell_metadata = cells, gene_metadata = gene_metadata)
#cid_cds <- estimate_size_factors(cid_cds)
#cid_cds <- classify_cells(cid_cds, classifier, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
##


## Automatec classification with garnett.
##
simplified_name <- gsub("#.*", "", names(classifier@signatures))
aggregated_markers <- list()
for (i in seq_along(simplified_name)){
    cell_type <- simplified_name[i]
    markers <- classifier@signatures[[i]]
    if (cell_type %in% names(aggregated_markers)){
        aggregated_markers[[cell_type]] <- c(aggregated_markers[[cell_type]], markers)
    } else {
        aggregated_markers[[cell_type]] <- markers
    }
}
aggregated_markers <- lapply(aggregated_markers, unique)
PanCancer_celltype <- data.frame(
    CellType = rep(names(aggregated_markers), sapply(aggregated_markers, length)),
    marker = unlist(aggregated_markers)
)
write.csv(PanCancer_celltype, "garnett_classifier/PanCancer_celltype_ref.csv", row.names = FALSE)
# generating a marker file.
output <- "garnett_classifier/garnett_marker_file.txt"
file_conn <- file(output, "w")
for (cell_type in names(aggregated_markers)){
    cell_type_format <- gsub(",", "", cell_type) #format cell type lines
    cat(paste0(">", cell_type_format, "\n"), file = file_conn)
    cat("expressed: ", paste(aggregated_markers[[cell_type]][1:9], collapse = ", "), "\n", file = file_conn)
}
close(file_conn)
# check marker
marker_file_path <- "garnett_classifier/garnett_marker_file.txt"
marker_check <- check_markers(cid_cds, marker_file_path, db = org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
garnett_markers <- plot_markers(marker_check)
ggsave("garnett_classifier/garnett_marker_plot.png", garnett_markers)
#saveRDS(cid_integrated, file = "cid_integrated.rds")
#cid_integrated <- readRDS("cid_integrated.rds")
# train classifier
set.seed(260)
filtered_checked_markers <- marker_check[marker_check$ambiguity <= 0.25, ]
filtered_checked_markers <- marker_check[marker_check$ambiguity <= 0.25, ]
PanCancer_classifier <- train_cell_classifier(cds = cid_cds,
    marker_file = "garnett_classifier/garnett_marker_file.txt",
    db = org.Hs.eg.db,
    cds_gene_id_type = "SYMBOL",
    num_unknown = 50,
    marker_file_gene_id_type = "SYMBOL"
)
feature_genes <- get_feature_genes(PanCancer_classifier, node = "root", db = org.Hs.eg.db, convert_ids = TRUE)
head(feature_genes)
## Now let's proceed to cell classification
pdata <- cid_integrated@meta.data
rownames(pdata) <- colnames(expm)
fdata <- data.frame(gene_short_name = rownames(expm))
rownames(fdata) <- rownames(expm)
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
cid_cds <- newCellDataSet(as(expm, "dgCMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
cid_cds <- estimateSizeFactors(cid_cds)
cid_cds_classified <- classify_cells(cid_cds, PanCancer_classifier, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
table(pData(cid_cds_classified)$cell_type)
table(pData(cid_cds_classified)$cluster_ext_type)
# all looking good, continue with clustring based on garnett classification
cid_cds_classified <- preprocess_cds(cid_cds_classified, num_dim = 150)







#spatial transcriptomics
dir = "filtered_feature_bc_matrix"
filter_matrix <- ReadMtx(
        mtx = paste0(dir, "/matrix.mtx.gz"), 
        features = paste0(dir, "/features.tsv.gz"), 
        cells = paste0(dir, "/barcodes.tsv.gz"), 
        feature.column = 1)
write10xCounts("filtered_feature_bc_matrix.h5", filter_matrix, type = "HDF5",
               genome = "hg38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))
CID44971.obj <- Load10X_Spatial(data.dir = "filtered_feature_bc_matrix")
Assays(CID44971.obj)
vln.plot <- VlnPlot(CID44971.obj, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(CID44971.obj , features = "nCount_Spatial") + theme(legend.position = "right")
vln.plot | count.plot



###genome = "hg38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))
CID44971.obj <- Load10X_Spatial(data.dir = "filtered_feature_bc_matrix")
Assays(CID44971.obj)
vln.plot <- VlnPlot(CID44971.obj, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(CID44971.obj , features = "nCount_Spatial") + theme(legend.position = "right")
vln.plot | count.plot
















genome = "hg38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))
CID44971.obj <- Load10X_Spatial(data.dir = "filtered_feature_bc_matrix")
Assays(CID44971.obj)
vln.plot <- VlnPlot(CID44971.obj, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(CID44971.obj , features = "nCount_Spatial") + theme(legend.position = "right")
vln.plot | count.plot














#integrate process and scale data individually.
sample_dir <- list.files(path = "data", pattern = "CID", full.names = TRUE)
sample_dir
cid <- list()
#Read the samples into the seurat object in a loop function
for (dir in sample_dir){
    project <- basename(dir)
    cid_matrix <- ReadMtx(
        mtx = paste0(dir, "/matrix.mtx.gz"), 
        features = paste0(dir, "/features.tsv.gz"), 
        cells = paste0(dir, "/barcodes.tsv.gz"), 
        feature.column = 1  
    )
    cid_matrix
    #cid_matrix = Read10X(data.dir = dir, gene.column = 1)


    cid_obj = CreateSeuratObject(counts = cid_matrix, min.cells = 3, min.features = 250)
    cid_obj[["percent.mt"]] <- PercentageFeatureSet(cid_obj, pattern = "^MT-")
    cid_obj <- subset(cid_obj, subset = nFeature_RNA > 200 & nCount_RNA > 250 & percent.mt < 20)
    cid_obj <- NormalizeData(cid_obj)
    cid_obj <- FindVariableFeatures(cid_obj, selection.method = "vst", nfeatures = 2000)
    cid_obj <- ScaleData(cid_obj)

    cid[[project]] <- cid_obj
}
vln_plot <- VlnPlot(cid_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "results/mt1_violin_plot.png", plot = vln_plot, width = 20, height = 6, dpi = 300)
#hvgs_all = SelectIntegrationFeatures(alldata.list)
cid_anchors <- FindIntegrationAnchors(object.list = cid, dims = 1:30)
cid_integrated <- IntegrateData(anchorset = cid_anchors, dims = 1:30) #end here, continue next line
#cid_integrated <- FindVariableFeatures(cid_integrated)
cid_integrated <- ScaleData(cid_integrated)
cid_integrated <- RunPCA(cid_integrated, npcs = 30)
DimPlot(cid_integrated, reduction = "pca") + NoLegend()

ElbowPlot(cid_integrated)
cid_integrated <- FindNeighbors(cid_integrated, dims = 1:30)
cid_integrated <- FindClusters(cid_integrated, resolution = 0.5)
cid_integrated <- RunUMAP(cid_integrated, dims = 1:30)
DimPlot(cid_integrated, reduction = "umap")
table(Idents(cid_integrated))
#markers_cluster0 <- FindMarkers(cid_integrated, ident.1 = 0)
#FindMarkers(cid_integrated)
head(Idents(cid_integrated), 5)
cluster5.markers <- FindMarkers(cid_integrated, ident.1 = 5, indent.2 = c(0, 3), min.pct = 0.25) #find markers for cluster 5 but very different from cluster 0 and 3
head(cluster5.marker, n = 5)
VlnPlot(cid_integrated, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2])) #plot the first two genes
all.markers <- FindAllMarkers(cid_integrated, anly.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x <- all.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(cid_integrated, features = x$gene[1:4])
FeaturePlot(cid_integrated, features = x$gene[5:8])
FeaturePlot(cid_integrated, features = x$gene[9:12])
FeaturePlot(cid_integrated, features = x$gene[13:15])
#plot the genes provided by the study
plot_markers <- FeaturePlot(cid_integrated, features = c("EPCAM", "MK167", "CD3D", "CD68", "MS4A1", "JCHAIN", "PECAM1", "PDGFRB"), combine = FALSE)
plot_markers <- lapply(X = plot_markers, FUN = function(x) x +
                                                            theme(plot.title = element_text(size = 8)) +
                                                            theme(axis.title.y = element_text(size = 5)) +
                                                            theme(axis.title.x = element_text(size = 5)) +
                                                            theme(axis.text.y = element_text(size =5)) +
                                                            theme(axis.text.x = element_text(size = 5)) +
                                                            theme(legend.positon = "none"))
CombinePlots(plots = plot_markers)

#annotate the clusters
new.cluster.ids <- c("Naive C04 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGRA+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(cid_integrated) 
cid_integrated <- RenameIdents(cid_integrated, new.cluster.ids)
DimPlot(cid_integrated, reduction = "pca", label = TRUE, pt.size = 0.9)
DimPlot(cid_integrated, reduction = "umap", label = TRUE)
levels(cid_integrated)
FeaturePlot(cid_integrated, features = c("EPCAM", "CD3D", "CD68", "MS4A1", 
                                         "JCHAIN", "PECAM1", "PDGFRB"))

new.cluster.ids <- c("Epithelial cells", "Proliferating cells", "T-cells", "Myeloid cells", 
                     "B-cells", "Plasmablasts", "Endothelial cells", "Mesenchymal cells", 
                     "Epithelial cells", "Proliferating cells", "T-cells", "Myeloid cells", 
                     "B-cells", "Plasmablasts", "Endothelial cells", "Mesenchymal cells")
names(new.cluster.ids) <- levels(cid_integrated)
cid_integratedd <- RenameIdents(cid_integrated, new.cluster.ids)
DimPlot(cid_integratedd, reduction = "pca", label = TRUE, pt.size = 0.9)
DimPlot(cid_integratedd, reduction = "umap", label = TRUE)











































#cid_integrated <- NormalizedData(cid_integrated)
#cid_integrated <- FindVariableFeatures(cid_integrated, selection.method ="vst", nfeatures = 2000)
#cid_integrated <- ScaleData(cid_integarated)
#DefaultAssay(cid_integrated) <- "integrated"

cid_integrated <- RunPCA(cid_integrated, features = VariableFeatures(object = cid_integrated), npcs = 30)
print(cid_integrated[["pca"]])
pca_plot <- DimPlot(cid_integrated, reduction = "pca") + NoLegend()
ggsave(filename = "PCA_plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)
pca_plot
#cid_integrated <- RunUMAP(cid_integrated, dims = 1:30)
cid_integrated <- RunPCA(cid_integrated, features = VariableFeatures(object = cid_integrated))






cell_ids <- names(cid)
merged_cid <- cid[[1]]
merged_cid <- merge(merged_cid, y = c(cid[2:10]), project = "CID")
merged_cid <- JoinLayers(merged_cid)
merged_cid
head(merged_cid@meta.data)
merged_cid@assays



# Perform SCTransform normalization and regress out mitochondrial percentage
options(future.globals.maxSize = 5000 * 1024^2)  # Set a higher limit 
merged_cid <- SCTransform(merged_cid, vars.to.regress = "percent.mt", variable.features.n = 1000)
top10 <- head(VariableFeatures(call), 10)
plot1 <- VariableFeaturePlot(call)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
# Find variable features



merged_cid[["percent.mt"]] <- PercentageFeatureSet(merged_cid, pattern = "^MT-")
vln_plot <- VlnPlot(merged_cid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "results/mt_violin_plot.png", plot = vln_plot, width = 20, height = 6, dpi = 300)

plot1 <- FeatureScatter(merged_cid, feature1 =  "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_cid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatter_plot <- plot1 + plot2
ggsave(filename = "results/feature_scatter_plot.png", plot = scatter_plot, width = 20, height = 6, dpi = 300)

merged_cid <- subset(merged_cid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt)
merged_cid <- FindVariableFeatures(merged_cid, selection.method = "vst", nfeatures = 2000)
merged_cid <- NormalizeData(merged_cid, normalization.method = "LogNormalize", scale.factor = 10000)
top10 <- head(VariableFeatures(merged_cid), 10)
plot1 <- VariableFeaturePlot(merged_cid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variable_feature_plot <- plot1 + plot2
ggsave(filename = "results/variable_feature_plot.png", plot = variable_feature_plot, width = 20, height = 6, dpi = 300)








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
