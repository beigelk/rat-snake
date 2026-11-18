###############################################################################
# Rule:     dotplot
# Purpose:  Dot plots of one marker set across one population's clusters at
#           one resolution
#
# Description:
#   The plotting code of the original dotplot.R is unchanged: a Seurat
#   DotPlot and a scCustomize Clustered_DotPlot, same sizes and colours.
#   Fixes: `assay_name`, `donefile` and `r_library` were never defined, the
#   rule wrote an empty .done file instead of the plot, and the clustered dot
#   plot failed when a marker set had fewer features than clusters. Features missing
#   from the assay are logged; a marker set with none left gets a placeholder
#   figure instead of failing the run.
#
# Inputs:   seurat_obj (compiled population object), features_file (filtered
#           marker table or gene list)
# Outputs:  png (DotPlot; its PDF and the DotPlotClustered_k<n> figures are
#           written next to it)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(qs2)
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(ComplexHeatmap)
  library(jcolors)
})
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "marker_features.R"))

p <- smk@params
out_base <- sub("\\.png$", "", smk@output$png)
out_dir <- dirname(out_base)

# -- FUNCTIONS -----------------------------------------------------------------
#' Seurat DotPlot and scCustomize Clustered_DotPlot of a marker set
#'
#' @param seu.obj Seurat object.
#' @param cluster.col Character. Cluster column to group by.
#' @param features Character vector of features in the default assay.
#' @param output.dir Character. Folder for the figures.
#' @param prefix Character. Start of the file names.
#'
#' @return NULL. Writes {prefix}__DotPlot and {prefix}__DotPlotClustered_k{n} (PDF and PNG).
plot_marker_dotplots <- function(seu.obj, cluster.col, features, output.dir, prefix) {
  Idents(seu.obj) <- cluster.col
  n_clusters <- length(unique(seu.obj@meta.data[[cluster.col]]))

  dp <- Seurat::DotPlot(seu.obj, group.by = cluster.col, features = features) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = sprintf("%s (%s)", prefix, cluster.col))

  plot_pdf_and_png(
    plot = dp,
    filename = file.path(output.dir, paste0(prefix, "__DotPlot")),
    width = 6 + n_clusters * 1.1,
    height = 4 + n_clusters * 0.1
  )

  # Clustered_DotPlot needs two or more clusters and features
  if (n_clusters < 2 || length(features) < 2) {
    log_message("Clustered dot plot skipped: needs two or more clusters and features", "WARN")
    return(invisible(NULL))
  }

  ht_opt$HEATMAP_LEGEND_PADDING <- unit(1, "cm")
  ht_opt$ANNOTATION_LEGEND_PADDING <- unit(1, "cm")
  ht_opt$DIMNAME_PADDING <- unit(2, "mm")
  ht_opt$heatmap_column_title_gp <- gpar(fontsize = 6, fontface = "plain")

  # k-means splits the features into k groups, so k must stay below the
  # number of features (the original always used k = clusters and failed
  # on short gene lists)
  k_val <- max(1, min(n_clusters, length(features) - 1))
  if (k_val < n_clusters) log_message(sprintf("Clustered dot plot: k = %d (fewer features than clusters)", k_val), "INFO")
  pal <- colorRampPalette(jcolors(palette = "rainbow"))(n_clusters)

  # The clustered plot is extra to the rule's output: if k-means cannot split
  # the features (e.g. all have the same expression), log it and go on
  cdp <- tryCatch(Clustered_DotPlot(
    seu.obj,
    group.by = cluster.col,
    features = features,
    cluster_ident = FALSE,
    show_ident_colors = TRUE,
    show_annotation_name = FALSE,
    legend_position = "right",
    row_names_side = "left",
    show_ident_legend = FALSE,
    colors_use_idents = pal,
    k = k_val,
    plot_km_elbow = FALSE,
    column_label_size = 10,
    row_label_size = 14,
    x_lab_rotate = 90,
    flip = TRUE
  ), error = function(e) {
    log_message(sprintf("Clustered dot plot skipped: %s", conditionMessage(e)), "WARN")
    NULL
  })
  if (is.null(cdp)) {
    ht_opt(RESET = TRUE)
    return(invisible(NULL))
  }

  plot_pdf_and_png(
    plot = cdp,
    filename = file.path(output.dir, paste0(prefix, "__DotPlotClustered_k", k_val)),
    width = 3 + n_clusters * 1.1,
    height = 3 + n_clusters * (if (n_clusters > 10) (0.1 + 1 / n_clusters) else 0.08)
  )

  ht_opt(RESET = TRUE)
  invisible(NULL)
}

# -- MAIN ----------------------------------------------------------------------
seu <- qs2::qs_read(smk@input$seurat_obj)
DefaultAssay(seu) <- p$assay

features <- present_features(seu, resolve_features(p$spec, smk@input$features_file), p$assay)
log_message(sprintf("%d features: %s", length(features), paste(features, collapse = ", ")), "INFO")

if (length(features) == 0) {
  placeholder_plot(out_base, p$prefix, "No features to plot (no markers passed the filter, or none are in the assay)")
} else {
  plot_marker_dotplots(seu, p$cluster_col, features, out_dir, p$prefix)
}
log_message("Dot plots written", "DONE")
