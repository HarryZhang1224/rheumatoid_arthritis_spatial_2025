figure_default_paths <- list(
  xenium_dir = "/data1/deyk/harry/RA_Xenium/data/post_qc_data_baysor",
  kmeans_result_file = "/data1/deyk/harry/RA_Xenium/results/nicheassay/kmeans_res/k_5_res.qs"
)

base_figure_packages <- c(
  "Seurat",
  "nicheDE",
  "qs",
  "glue",
  "ggplot2",
  "patchwork",
  "latex2exp",
  "scCustomize",
  "scico",
  "ComplexHeatmap",
  "ggsci",
  "ggrepel",
  "org.Hs.eg.db",
  "clusterProfiler",
  "RColorBrewer",
  "Polychrome",
  "stringr",
  "parallel",
  "ggpubr",
  "dplyr",
  "tidyr"
)

load_figure_libraries <- function(extra_packages = character()) {
  packages <- unique(c(base_figure_packages, extra_packages))
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

colors_all <- c(
  "#DA80DA", "#815481", "#C040C0", "#E1AFE1", "#3F0034",
  "#EDABB9", "#EB5C79", "#A06A75", "#C00028",
  "#EB675E", "#A23E36",
  "#540F54", "#53407F",
  "#DFA38A", "#8C3612", "#623623", "#916350", "#DAC3C3",
  "#F8770B", "#E09E3A", "#CD7225", "#FFC990", "#AC5812",
  "#FEE083", "#897538", "#E7B419", "#BCA048",
  "#6F8BE2", "#3053BC",
  "#6D9F58", "#9EB766", "#BDCB10", "#3A6527", "#9EA743",
  "#E2E8A7", "#5A6209", "#8FE36B",
  "#818A31",
  "#9FC5E8", "#23D9F1",
  "#64C6A6",
  "#AAAAAA"
)

ct_cols <- c(
  "M2" = colors_all[10],
  "M5" = colors_all[28],
  "F2" = colors_all[30],
  "F7" = colors_all[24],
  "M1" = colors_all[4],
  "NK" = colors_all[39],
  "M10" = colors_all[17],
  "M3" = colors_all[32],
  "F3" = colors_all[29],
  "F4" = colors_all[12],
  "B3" = colors_all[41],
  "E0" = colors_all[19],
  "B1" = colors_all[7],
  "B2" = colors_all[11]
)

height_width_param <- c("6&8", "6.5&7", "6&8", "5&9", "6&8", "7.5&4", "6&10", "6&10", "8&10", "6&8")
names(height_width_param) <- c("RA401", "RA331", "RA442", "RA443", "RA457", "RA480", "RA494", "RA519", "RA362", "RA489")

get_all_sample_ids <- function(data_dir = figure_default_paths$xenium_dir) {
  all_ra_samples <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
  stringr::str_extract(all_ra_samples, "RA\\d+[A-Z]*")
}

get_height_width_param <- function() {
  height_width_param
}

get_sample_dimensions <- function(sample_id, dimensions = get_height_width_param()) {
  as.numeric(stringr::str_split(dimensions[[sample_id]], "&")[[1]])
}

load_xenium_sample <- function(sample_id, data_dir = figure_default_paths$xenium_dir) {
  qs::qread(glue::glue("{data_dir}/{sample_id}/final_version_seurat_obj.qs"))
}

load_kmeans_niche_df <- function(kmeans_path = figure_default_paths$kmeans_result_file) {
  kmeans_res <- qs::qread(kmeans_path)$cluster
  kmeans_niche_df <- data.frame(Niche = kmeans_res) %>%
    dplyr::mutate(sample = gsub("_.+", "", names(kmeans_res)))
  kmeans_niche_df <- split(kmeans_niche_df, kmeans_niche_df$sample)
  lapply(kmeans_niche_df, function(df) {
    rownames(df) <- gsub("RA[0-9]+_", "", rownames(df))
    df %>%
      dplyr::mutate(Niche = factor(Niche))
  })
}

prepare_xenium_obj_for_merge <- function(obj, add_coords = TRUE, clear_images = TRUE) {
  DefaultAssay(obj) <- "Xenium"
  for (assay_name in c("SCT", "ControlCodeword", "ControlProbe")) {
    if (assay_name %in% names(obj@assays)) {
      obj[[assay_name]] <- NULL
    }
  }
  if (add_coords) {
    coords <- GetTissueCoordinates(obj)
    rownames(coords) <- coords$cell
    coords <- dplyr::select(coords, -cell)
    obj <- AddMetaData(obj, coords)
  }
  if (clear_images) {
    obj@images <- list()
  }
  obj
}

load_prepared_xenium_samples <- function(sample_ids = get_all_sample_ids(),
                                         data_dir = figure_default_paths$xenium_dir,
                                         attach_niches = FALSE,
                                         kmeans_niche_df = NULL,
                                         add_coords = TRUE,
                                         clear_images = TRUE) {
  all_xenium_samples <- setNames(lapply(sample_ids, load_xenium_sample, data_dir = data_dir), sample_ids)
  all_xenium_samples <- lapply(all_xenium_samples, prepare_xenium_obj_for_merge, add_coords = add_coords, clear_images = clear_images)

  if (attach_niches) {
    if (is.null(kmeans_niche_df)) {
      kmeans_niche_df <- load_kmeans_niche_df()
    }
    all_xenium_samples <- lapply(names(all_xenium_samples), function(sample_id) {
      AddMetaData(all_xenium_samples[[sample_id]], kmeans_niche_df[[sample_id]])
    })
    names(all_xenium_samples) <- sample_ids
  }

  all_xenium_samples
}

load_merged_xenium_samples <- function(sample_ids = get_all_sample_ids(),
                                       data_dir = figure_default_paths$xenium_dir,
                                       attach_niches = TRUE,
                                       kmeans_niche_df = NULL,
                                       add_coords = TRUE,
                                       clear_images = TRUE,
                                       normalize = TRUE) {
  if (attach_niches && is.null(kmeans_niche_df)) {
    kmeans_niche_df <- load_kmeans_niche_df()
  }

  all_xenium_samples <- load_prepared_xenium_samples(
    sample_ids = sample_ids,
    data_dir = data_dir,
    attach_niches = attach_niches,
    kmeans_niche_df = kmeans_niche_df,
    add_coords = add_coords,
    clear_images = clear_images
  )

  if (length(all_xenium_samples) == 1) {
    all_xenium_merged <- all_xenium_samples[[1]]
  } else {
    all_xenium_merged <- merge(all_xenium_samples[[1]], all_xenium_samples[-1])
  }

  all_xenium_merged[["Xenium"]] <- JoinLayers(all_xenium_merged[["Xenium"]])
  if (normalize) {
    all_xenium_merged <- NormalizeData(all_xenium_merged, normalization.method = "LogNormalize")
  }

  all_xenium_merged
}
