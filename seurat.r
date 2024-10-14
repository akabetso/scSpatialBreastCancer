
library(Seurat)


sample_dir <- list.files(path = "data", pattern = "CID", full.names = TRUE)

seurat_objects <- list()
#Read the samples into the seurat object in a loop function
for (dir in sample_dir){
    #Read te count matrix from the directory
    count_matrix = Read10X(data.dir = dir)
    #create seurat object
    seurat_obj <- CreateSeuratObject(counts = count_matrix, project = dir)
    #append the object to the empty list
    seurat_objects[[dir]] <- seruat_obj

}

sc_data = Read10X(data.dir = "data/CID3586")

ls data/CID3586
