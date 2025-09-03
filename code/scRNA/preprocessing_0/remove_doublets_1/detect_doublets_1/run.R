library(Seurat)
library(DoubletFinder)
library(qs)
library(glue)

sample_id <- commandArgs(trailingOnly = TRUE)
doubletfinder_intermediate <- qread(glue("/data1/deyk/harry/RA_Xenium/results/DoubletFinder/{sample_id}.qs"))
optim_pk <- as.numeric(as.character(doubletfinder_intermediate$bcmvn[, "pK"][which.max(doubletfinder_intermediate$bcmvn[, "BCmetric"])]))
seurat_obj <- qread(glue("/data1/deyk/harry/RA_Xenium/data/Ian_old_scRNA_07062025/seurat/{sample_id}.qs"))
seurat_obj <- FindNeighbors(seurat_obj, dims=1:30, verbose=FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution=0.2)
homotypic_proportion <- modelHomotypic(seurat_obj$seurat_clusters) # homotypic proportion estimate
nExp_poi <- round(0.05*ncol(seurat_obj)) # doublet formation rate is 0.05 according to 10X when loading 10k cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic_proportion)) # heterotypic doublet rate
message(glue("Optimal pK: {optim_pk}, expected number of heterotypic doublets: {nExp_poi.adj}"))
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = optim_pk, nExp = nExp_poi.adj, sct = TRUE)
qsave(seurat_obj, glue("/data1/deyk/harry/RA_Xenium/data/Ian_old_scRNA_07062025/doublet_removed_seurat/{sample_id}.qs"))