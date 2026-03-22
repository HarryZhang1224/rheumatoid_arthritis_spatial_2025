mid_rescaler <- function(mid) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}

spatialVariablePlot <- function(obj, plot_cts, focal_ct = NULL, focal_size = 0.85, color_map = ct_cols) {
  coords <- GetTissueCoordinates(obj)
  rownames(coords) <- coords$cell
  ap_ratio <- diff(range(coords[, "x"])) / diff(range(coords[, "y"]))
  scale_bar_x_start <- min(coords[, "y"]) + 50
  scale_bar_x_end <- min(coords[, "y"]) + 50 + 1000
  scale_bar_y_pos <- min(coords[, "x"] + 50)
  coords <- dplyr::select(coords, -cell)
  cells_of_interest <- which(Idents(obj) %in% plot_cts)
  coords$celltype <- "Background"
  coords$celltype[cells_of_interest] <- as.character(Idents(obj)[cells_of_interest])
  coords$celltype <- factor(coords$celltype, levels = c("Background", sort(unique(as.character(Idents(obj)[cells_of_interest])))))
  plot_order <- seq_len(nlevels(coords$celltype))
  names(plot_order) <- levels(coords$celltype)
  coords$plot_order <- plot_order[as.character(coords$celltype)]
  coords <- dplyr::arrange(coords, plot_order)

  alpha_vals <- rep(1, nlevels(coords$celltype))
  names(alpha_vals) <- levels(coords$celltype)
  alpha_vals["Background"] <- 0.2

  size_vals <- rep(0.4, nlevels(coords$celltype))
  names(size_vals) <- levels(coords$celltype)
  if (length(focal_ct) > 0) {
    size_vals[focal_ct] <- focal_size
  }

  ggplot(coords, aes(y, x, fill = celltype, alpha = celltype, size = celltype)) +
    geom_point(stroke = 0.05, shape = 21, color = "black") +
    scale_fill_manual(values = color_map, breaks = setdiff(levels(coords$celltype), "Background")) +
    scale_alpha_manual(values = alpha_vals, breaks = setdiff(levels(coords$celltype), "Background"), name = "") +
    scale_size_manual(values = size_vals, breaks = setdiff(levels(coords$celltype), "Background"), name = "") +
    theme_void() +
    theme(
      legend.key.size = unit(0.3, "cm"),
      aspect.ratio = ap_ratio,
      legend.text = element_text(size = 4),
      plot.title = element_text(size = 6)
    ) +
    labs(color = "", fill = "", size = "") +
    guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    annotate("segment", x = scale_bar_x_start, xend = scale_bar_x_end, y = scale_bar_y_pos, yend = scale_bar_y_pos, linewidth = 0.5) +
    annotate("text", x = (scale_bar_x_start + scale_bar_x_end) / 2, y = scale_bar_y_pos, label = paste0(1, " mm"), vjust = 1.5, size = 1.5)
}

spatialClusterPlotAllCT <- function(df,
                                    ct_column,
                                    x_column = "centroid_x",
                                    y_column = "centroid_y",
                                    show.legend = TRUE,
                                    stroke = 0.1,
                                    size = 0.6,
                                    palette = paletteer::paletteer_d("pals::glasbey"),
                                    title = "") {
  scale_bar_x_start <- min(df[, x_column, drop = TRUE]) + 50
  scale_bar_x_end <- min(df[, x_column, drop = TRUE]) + 50 + 1000
  scale_bar_y_pos <- min(df[, y_column, drop = TRUE] + 50)

  ggplot(df, aes(x = !!sym(x_column), y = !!sym(y_column), fill = !!sym(ct_column))) +
    geom_point(stroke = stroke, size = size, shape = 21, show.legend = show.legend) +
    theme_void() +
    theme(
      legend.key.size = unit(0.2, "cm"),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, face = "bold")
    ) +
    coord_fixed() +
    scale_fill_manual(values = palette, na.translate = FALSE) +
    guides(fill = guide_legend(override.aes = list(size = 1.5))) +
    labs(x = "sp1", y = "sp2", fill = "") +
    annotate("segment", x = scale_bar_x_start, xend = scale_bar_x_end, y = scale_bar_y_pos, yend = scale_bar_y_pos, linewidth = 0.5) +
    annotate("text", x = (scale_bar_x_start + scale_bar_x_end) / 2, y = scale_bar_y_pos, label = paste0(1, " mm"), vjust = 1.5, size = 1.5) +
    ggtitle(title)
}

plotUMAPFeature <- function(df,
                            embedding_1,
                            embedding_2,
                            variable,
                            palette = rev(paletteer::paletteer_c("grDevices::Reds", n = 50)),
                            stroke = 0.03,
                            cap = FALSE,
                            cap_quantile = 0.98,
                            size = 0.3,
                            order_variable = variable,
                            legend_title = variable,
                            midpoint = NA,
                            title = "",
                            order = FALSE) {
  rescaler <- if (is.na(midpoint)) scales::rescale else mid_rescaler(midpoint)
  if (cap) {
    df[[variable]] <- pmin(df[[variable]], stats::quantile(df[[variable]], cap_quantile, na.rm = TRUE))
  }
  if (order) {
    df <- df %>% dplyr::arrange(!!sym(order_variable))
  }

  ggplot(df, aes(!!sym(embedding_1), !!sym(embedding_2), fill = !!sym(variable))) +
    geom_point(size = size, stroke = stroke, shape = 21) +
    scale_fill_gradientn(colors = palette, rescale = rescaler) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 8),
      legend.key.size = unit(3, "mm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    ) +
    ggtitle(title) +
    labs(fill = legend_title)
}

spatialValuePlot <- function(df,
                             variable,
                             x_column = "centroid_x",
                             y_column = "centroid_y",
                             pt.size = 0.2,
                             stroke = 0.05,
                             ct.column = "",
                             midpoint = NA,
                             plt.ct = "",
                             title = variable,
                             palette = paletteer::paletteer_c("viridis::inferno", n = 30),
                             order = TRUE,
                             legend_title = variable,
                             custom_scale = FALSE,
                             scale_externel = NULL) {
  rescaler <- if (is.na(midpoint)) scales::rescale else mid_rescaler(midpoint)
  sp_aspect_ratio <- diff(range(df[, y_column])) / diff(range(df[, x_column]))
  if (nchar(plt.ct) > 0 && nchar(ct.column) > 0) {
    df <- df %>%
      dplyr::filter(!!sym(ct.column) == plt.ct)
  }
  if (order) {
    df <- df %>%
      dplyr::arrange(!!sym(variable))
  }
  color_bar <- if (custom_scale) scale_externel else scale_fill_gradientn(colors = palette, rescale = rescaler)
  scale_bar_x_start <- min(df[, x_column, drop = TRUE]) + 50
  scale_bar_x_end <- min(df[, x_column, drop = TRUE]) + 50 + 1000
  scale_bar_y_pos <- min(df[, y_column, drop = TRUE] + 50)

  ggplot(df, aes(!!sym(x_column), !!sym(y_column), fill = !!sym(variable))) +
    geom_point(stroke = stroke, size = pt.size, shape = 21) +
    theme_bw() +
    theme(
      aspect.ratio = sp_aspect_ratio,
      legend.key.size = unit(0.2, "cm"),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 5),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, face = "bold")
    ) +
    color_bar +
    labs(x = "sp1", y = "sp2", fill = legend_title) +
    ggtitle(title) +
    annotate("segment", x = scale_bar_x_start, xend = scale_bar_x_end, y = scale_bar_y_pos, yend = scale_bar_y_pos, linewidth = 0.5) +
    annotate("text", x = (scale_bar_x_start + scale_bar_x_end) / 2, y = scale_bar_y_pos, label = paste0(1, " mm"), vjust = 1.5, size = 1.5)
}

getSharedScales <- function(values,
                            midpoint,
                            palette = rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging", n = 100))) {
  scale_fill_gradientn(
    colors = palette,
    limits = range(values, na.rm = TRUE),
    rescaler = if (is.na(midpoint)) scales::rescale else mid_rescaler(midpoint)
  )
}
