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
library(SpatialExperiment)
library(dplyr)
library(tidyr)

input_args <- commandArgs(trailingOnly = TRUE)
sample_id <- input_args[1]
radius <- input_args[2]
obj <- qread(glue("/data1/deyk/harry/RA_Xenium/data/post_qc_data_baysor/{sample_id}/final_version_seurat_obj.qs"))
obj@meta.data$celltype_subcluster <- factor(obj@meta.data$celltype_subcluster)
nn_mat <- qread(glue("/data1/deyk/harry/RA_Xenium/data/neighbors/neighbor_matrices/Baysor/{sample_id}_neighbors_radius_{radius}_final.qs"))

ct_counts <- table(obj@meta.data$celltype_subcluster)
# Only look at cells with more than 30 counts in the dataset, otherwise the estimate is too noisy
all_unique_cell_types <- names(ct_counts)[ct_counts > 30]

# Get the number of neighbors for each cell
num_neighbors <- apply(nn_mat, 1, function(neighbors){
    sum(neighbors > 0) - 1
})

# Only sample cells with neighbors
cells_to_sample_from <- which(num_neighbors > 0)

# Get the number of cells to sample for each cell type
sampling_sizes <- sapply(ct_counts[all_unique_cell_types], function(count){
    min(count, length(cells_to_sample_from))
})

# Get the cell type assignment vector
ct_assignment_vector <- (obj@meta.data$celltype_subcluster)

CLQ_bg_sim_all_cts <- mclapply(all_unique_cell_types, function(target_ct){
    # Pull 100 CLQ scores under the null assumption, returns a cell type X number of iterations CLQ value matrix
    CLQ_bg_sim <- do.call(cbind, lapply(seq_len(100), function(iter){
        # At each iteration (seed)
        set.seed(iter)
        # Assign target cell types to random locations in the tissue
        sampled_cells <- sample(cells_to_sample_from, size=sampling_sizes[target_ct], replace=FALSE)
        # Returns a cell type X number of sampled cells Cba values
        Cba_bg <- do.call(cbind, lapply(sampled_cells, function(idx){
            # For each sampled cell
            # Extract the neighbors
            community <- nn_mat[idx, ]
            neighbor_cell_idx <- community[community > 0]
            # Remove target cell (self) from the neighborhood
            neighbor_cell_idx <- neighbor_cell_idx[-1]
            # Total number of cells in the neighborhood
            denom <- max(length(neighbor_cell_idx), 1)
            # Get the cell types of the neighbors and their proportions in the neighborhood of the target cell
            table(ct_assignment_vector[neighbor_cell_idx])/denom
        }))
        # Get the CLQ value for each cell type = observed proportions of infiltrating cell types in target cell type neighborhood / expected proprotion of infiltrating cell types based on tissue composition
        CLQ_bg <- (rowSums(Cba_bg)/sampling_sizes[target_ct])/(ct_counts[rownames(Cba_bg)]/(ncol(obj) - 1))
        CLQ_bg
    }))
    return(CLQ_bg_sim)
}, mc.cores=35)
names(CLQ_bg_sim_all_cts) <- all_unique_cell_types

qsave(CLQ_bg_sim_all_cts, glue("/data1/deyk/harry/RA_Xenium/results/CLQ_res/Baysor/null_simulated_CLQs/{sample_id}_radius_{radius}_simulated_CLQs.qs"))
