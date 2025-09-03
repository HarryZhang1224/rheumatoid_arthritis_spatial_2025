library(Seurat)
library(glue)
library(qs)
library(nicheDE)
library(dplyr)
library(tidyr)

source("/data1/deyk/harry/RA_Xenium/job_scripts/runNicheDE/NicheDECustomGamma.R")
sample_id <- commandArgs(trailingOnly = TRUE)
obj <- qread(glue("/data1/deyk/harry/RA_Xenium/data/post_qc_data_baysor/{sample_id}/final_version_seurat_obj.qs"))
obj@meta.data$celltype_subcluster <- as.character(obj@meta.data$celltype_subcluster)
obj@meta.data$celltype_subcluster[grepl("^M[25]", obj@meta.data$celltype_subcluster)] <- "M25"
obj@meta.data$celltype_subcluster <- factor(obj@meta.data$celltype_subcluster)
Idents(obj) <- "celltype_subcluster"
ct_counts <- table(Idents(obj))
cts_to_test <- names(ct_counts)[ct_counts > 50]
cells_to_keep <- which(Idents(obj) %in% cts_to_test)
obj <- subset(obj, cells=colnames(obj)[cells_to_keep])
raw_count_mat <- t(obj[["Xenium"]]$counts)
ct_mat <- obj@meta.data %>% 
    mutate(cell=rownames(.), celltype_subcluster=as.character(celltype_subcluster)) %>%
    select(cell, celltype_subcluster)
coords <- GetTissueCoordinates(obj)
rownames(coords) <- coords$cell
coords <- coords %>% select(-cell)
lib_mat <- CreateLibraryMatrix(raw_count_mat, cell_type=ct_mat)
ct_mat %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = celltype_subcluster, values_from = value, values_fill = list(value = 0)) %>%
    data.frame(., check.names=FALSE) -> ct_mat
rownames(ct_mat) <- ct_mat$cell
ct_mat %>% select(-cell) -> ct_mat
ct_mat <- as.matrix(ct_mat[, rownames(lib_mat)])
# The sigma values have been checked to include respectively on average 9, 20, and 30 neighbors
NDE_obj = CreateNicheDEObject(raw_count_mat,
                              coords,
                              lib_mat,
                              ct_mat,
                              sigma=c(175,250,300))
NDE_obj = CalculateEffectiveNicheLargeScale(NDE_obj,batch_size = 1000, cutoff = 0.05)
gamma_per_ct <- rep(0.8, nrow(NDE_obj@ref_expr))
names(gamma_per_ct) <- rownames(NDE_obj@ref_expr)
gamma_per_ct[grepl("^[FM]", names(gamma_per_ct))] <- 0.7
NDE_obj = niche_DE_no_parallel(NDE_obj,
                   C=150,
                   M=50,
                   gamma=gamma_per_ct,
                   print=T,
                   Int=T,
                   batch=F,
                   self_EN=T)
qsave(NDE_obj, glue("/data1/deyk/harry/RA_Xenium/results/nicheDE/{sample_id}_nicheDE_all_subtype_M25_combined.qs"))                           