library(Seurat,lib.loc="/lila/data/deyk/harry/spatial/conda_envs/R_updated/lib/R/seurat_10X_100324")
library(arrow)
library(qs)
library(reticulate)
library(dplyr)
library(glue)

use_condaenv("ST_env", required = TRUE)
options(future.globals.maxSize = 8.0 * 1e9)
sample_id <- commandArgs(trailingOnly = TRUE)
data_dir <- glue("/data/deyk/harry/spatial/RA_Xenium/data/xenium_ranger_baysor/{sample_id}/outs/")
xenium_obj <- LoadXenium(data_dir, fov="fov")
median_nfeature <- median(xenium_obj@meta.data$nFeature_Xenium)
mad_nfeature <- mad(xenium_obj@meta.data$nFeature_Xenium)
median_ncount <- median(xenium_obj@meta.data$nCount_Xenium)
mad_ncount <- mad(xenium_obj@meta.data$nCount_Xenium)
nfeature_upper_lim <- median_nfeature + 3*mad_nfeature
nfeature_lower_lim <- 10
ncount_upper_lim <- median_ncount + 3*mad_ncount
ncount_lower_lim <- 15
# Filter cells based on both sample-specific and global thresholds
xenium_obj <- subset(xenium_obj, subset=nFeature_Xenium >= nfeature_lower_lim & nFeature_Xenium < nfeature_upper_lim & nCount_Xenium >= ncount_lower_lim & nCount_Xenium < ncount_upper_lim)
xenium_obj <- SCTransform(xenium_obj, assay = "Xenium")
xenium_obj <- RunPCA(xenium_obj, npcs=50)
pct <- xenium_obj[["pca"]]@stdev / sum(xenium_obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1] # This finds the elbow
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
pcs <- min(co1, co2)

message(glue("Number of PCs: {pcs}"))

xenium_obj <- RunUMAP(xenium_obj, dims=1:pcs)
xenium_obj <- FindNeighbors(xenium_obj, reduction = "pca", dims = 1:pcs)
xenium_obj <- FindClusters(xenium_obj, resolution = 0.8, algorithm=4, method="igraph")

system(glue("mkdir -p /lila/data/deyk/harry/spatial/RA_Xenium/data/post_qc_data_baysor/{sample_id}"))

qsave(xenium_obj, glue("/lila/data/deyk/harry/spatial/RA_Xenium/data/post_qc_data_baysor/{sample_id}/QCed_seurat_obj.qs"))