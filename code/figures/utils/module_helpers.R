getWeightedModuleScore <- function(features,
                                   feature_weights,
                                   seurat_obj,
                                   module_name,
                                   assay = "RNA",
                                   layer = "data",
                                   ctrl_size = 100) {
  assay.data <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  data.avg <- rowMeans(assay.data)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(
    x = data.avg + stats::rnorm(n = length(data.avg)) / 1e+30,
    n = 24,
    labels = FALSE,
    right = FALSE
  )
  names(data.cut) <- names(data.avg)
  names(feature_weights) <- features
  feature_gene_bins <- data.cut[features]
  ctrl.genes <- unlist(lapply(feature_gene_bins, function(bin) {
    names(sample(data.cut[data.cut == bin], size = ctrl_size, replace = FALSE))
  }))
  ctrl.genes.weights <- rep(feature_weights, each = ctrl_size)
  dup_tbl <- tapply(ctrl.genes.weights, ctrl.genes, sum)
  ctrl.genes <- names(dup_tbl)
  ctrl.genes.weights <- as.numeric(dup_tbl)
  feature_scores <- as.numeric(Matrix::crossprod(assay.data[features, ], matrix(feature_weights, nrow = length(feature_weights)))) / sum(abs(feature_weights))
  ctrl_scores <- as.numeric(Matrix::crossprod(assay.data[ctrl.genes, ], matrix(ctrl.genes.weights, nrow = length(ctrl.genes.weights)))) / sum(abs(ctrl.genes.weights))
  module_score <- feature_scores - ctrl_scores
  module_score <- data.frame(module_score)
  colnames(module_score) <- module_name
  rownames(module_score) <- colnames(assay.data)
  module_score
}

getWeightedModuleScorescRNA <- function(features,
                                        feature_weights,
                                        seurat_obj,
                                        module_name) {
  getWeightedModuleScore(
    features = features,
    feature_weights = feature_weights,
    seurat_obj = seurat_obj,
    module_name = module_name,
    assay = "RNA",
    layer = "data",
    ctrl_size = 100
  )
}

organizeDEGsList <- function(fc_cutoff,
                             pval_cutoff,
                             degs_list,
                             sc_ref) {
  if (fc_cutoff < 0) {
    degs.tmp <- degs_list %>%
      dplyr::filter(logFC <= fc_cutoff & P.Value <= pval_cutoff)
  } else {
    degs.tmp <- degs_list %>%
      dplyr::filter(logFC >= fc_cutoff & P.Value <= pval_cutoff)
  }

  degs.tmp <- degs.tmp %>%
    dplyr::mutate(gene = rownames(.)) %>%
    dplyr::filter(gene %in% rownames(sc_ref)) %>%
    dplyr::arrange(dplyr::desc(abs(logFC)))

  setNames(abs(degs.tmp$logFC), degs.tmp$gene)
}
