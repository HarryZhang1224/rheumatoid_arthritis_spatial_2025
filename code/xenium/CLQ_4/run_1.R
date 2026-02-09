library(Seurat)
library(glue)
library(qs)
library(nicheDE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(parallel)
library(ggplot2)
library(stringr)
library(ggsci)
library(reshape2)
library(scales)
library(hues)
library(RColorBrewer)
library(latex2exp)
library(MASS)
library(irr)
library(patchwork)
library(ggdark)
library(CellChat)
library(spdep)
library(dplyr)
library(tidyr)

# computeCLQ <- function(cell_metadata, neighbors, Ca, Cb){
#     message(glue("Computing CLQ for target cell type {Ca} and infiltrating cell type {Cb}"))
#     # cell_metadata: metadata that has information about each cell's cluster id
#     # neighbors: matrix that has the cell ids within the neighborhood for each cell
#     # Ca: target cell type name
#     # Cb: infiltrating cell type name
#     num_neighbors <- apply(neighbors, 1, function(neighbors){
#             sum(neighbors > 0) - 1
#     })
#     infiltrating_cells <- which(cell_metadata$celltype_subcluster==Cb)
#     target_cells <- which(cell_metadata$celltype_subcluster==Ca)
#     Na <- length(target_cells)
#     Nb <- length(infiltrating_cells)
#     N <- nrow(cell_metadata)
#     # Compute Cba for each cell
#     Cba <- apply(matrix(neighbors[target_cells, ], ncol=ncol(neighbors)), 1, function(community){
#             denom <- max((sum(community > 0) - 1), 1)
#             length(intersect(community, infiltrating_cells))/denom
#     })
#     # This gives the total number of infiltrating cells in the neighborhood of all target cells
#     Cba_sum <- sum(Cba)
#     # Global CLQ value
#     CLQ_global <- (Cba_sum/Na)/(Nb/(N-1))
    
#     # Compute empirical p-value for the global CLQ by assigning target cells to random coordinates
#     cells_to_sampled_from <- which(num_neighbors > 0)
#     sampling_size <- min(length(target_cells), length(cells_to_sampled_from))
#     # Compute the global CLQ under the null
#     CLQ_global_background <- lapply(seq_len(100), function(iter){
#         set.seed(iter)
#         sampled_cells <- sample(cells_to_sampled_from, size=sampling_size, replace=FALSE)
#         Cba_global_background <- apply(matrix(neighbors[sampled_cells, ], ncol=ncol(neighbors)), 1, function(community){
#             denom <- max((sum(community > 0) - 1), 1)
#             length(intersect(community, infiltrating_cells))/denom
#         })
#         Cba_sum_global_background <- sum(Cba_global_background)
#         (Cba_sum_global_background/sampling_size)/(Nb/(N-1))
#     })
#     CLQ_global_background <- unlist(CLQ_global_background)
#     # Compute one-sided p-values
#     pval_global_CLQ <- mean(CLQ_global_background >= CLQ_global)
#     pval_global_CLQ_separation <- mean(CLQ_global_background <= CLQ_global)
    
#     return(list(CLQ_global=CLQ_global,
#                 CLQ_background=CLQ_global_background,
#                 pval_global_CLQ=pval_global_CLQ,
#                 pval_global_CLQ_separation=pval_global_CLQ_separation))
# }
computeCLQ <- function(cell_metadata, neighbors, Ca, Cb, null_bg){
    message(glue("Computing CLQ for target cell type {Ca} and infiltrating cell type {Cb}"))
    # cell_metadata: metadata that has information about each cell's cluster id
    # neighbors: matrix that has the cell ids within the neighborhood for each cell
    # Ca: target cell type name
    # Cb: infiltrating cell type name
    # null_bg: CLQ values under the null assumption
    num_neighbors <- apply(neighbors, 1, function(neighbors){
            sum(neighbors > 0) - 1
    })
    infiltrating_cells <- which(cell_metadata$celltype_subcluster==Cb)
    target_cells <- which(cell_metadata$celltype_subcluster==Ca)
    Na <- length(target_cells)
    Nb <- length(infiltrating_cells)
    N <- nrow(cell_metadata)
    # Compute Cba for each cell
    Cba <- apply(matrix(neighbors[target_cells, ], ncol=ncol(neighbors)), 1, function(community){
            denom <- max((sum(community > 0) - 1), 1)
            length(intersect(community, infiltrating_cells))/denom
    })
    # This gives the total number of infiltrating cells in the neighborhood of all target cells
    Cba_sum <- sum(Cba)
    # Global CLQ value
    CLQ_global <- (Cba_sum/Na)/(Nb/(N-1))
    
    # Compute empirical p-value by comparing to CLQ values under the null that cells are randomly distributed regardless of their identities
    CLQ_global_background <- null_bg[[Ca]][Cb, ]
    pval_global_CLQ <- mean(CLQ_global_background >= CLQ_global)
    pval_global_CLQ_separation <- mean(CLQ_global_background <= CLQ_global)
    
    return(list(CLQ_global=CLQ_global,
                CLQ_background=CLQ_global_background,
                pval_global_CLQ=pval_global_CLQ,
                pval_global_CLQ_separation=pval_global_CLQ_separation))
}
input_args <- commandArgs(trailingOnly = TRUE)
sample_id <- input_args[1]
radius <- input_args[2]
obj <- qread(glue("/data1/deyk/harry/RA_Xenium/data/post_qc_data_baysor/{sample_id}/final_version_seurat_obj.qs"))
nn_mat <- qread(glue("/data1/deyk/harry/RA_Xenium/data/neighbors/neighbor_matrices/Baysor/{sample_id}_neighbors_radius_{radius}_final.qs"))
ct_counts <- table(obj@meta.data$celltype_subcluster)
all_unique_cell_types <- names(ct_counts)[ct_counts > 30]
ct_pair_test_mat <- expand.grid(all_unique_cell_types, all_unique_cell_types, stringsAsFactors=FALSE)
null_clqs <- qread(glue("/data1/deyk/harry/RA_Xenium/results/CLQ_res/Baysor/null_simulated_CLQs/{sample_id}_radius_{radius}_simulated_CLQs.qs"))
clq_res_all <- mclapply(seq_len(nrow(ct_pair_test_mat)), function(idx){
    ca.tmp <-  ct_pair_test_mat$Var1[idx]
    cb.tmp <- ct_pair_test_mat$Var2[idx]
    computeCLQ(obj@meta.data, nn_mat, ca.tmp, cb.tmp, null_clqs)
}, mc.cores=30)
names(clq_res_all) <- paste(ct_pair_test_mat$Var1, ct_pair_test_mat$Var2, sep=",")
qsave(clq_res_all, glue("/data1/deyk/harry/RA_Xenium/results/CLQ_res/Baysor/{sample_id}_CLQ_final_{radius}.qs"))