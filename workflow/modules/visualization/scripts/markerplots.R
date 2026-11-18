###############################################################################
# Rule:     markerplots
# Purpose:  One combined figure per marker set and resolution: stacked
#           violins and a labelled cluster UMAP above feature UMAPs
#
# Description:
#   Wired from combo_plot.R, which looped over hard-coded marker folders and a
#   hard-coded population, and wrote to a 'V5' folder (and, for the combo
#   plot, to the working directory). Each job now plots one population,
#   resolution, assay and marker set from the config. The three plots and
#   their palettes, sizes and layout are unchanged:
#     __FeatureDim    feature UMAPs (scplotter FeatureStatPlot, "dim")
#     __StackViolin   stacked violins per cluster
#     __ComboPlot     (violins | cluster UMAP) / feature UMAPs
#   Fixes: `row2_t_margin` was used before it was set; cluster labels became
#   NA when the cluster column was not a factor; the UMAP name was fixed
#   (now visualization.reduction). A marker set with no features gets a
#   placeholder figure.
#
# Inputs:   seurat_obj (compiled population object), features_file (filtered
#           marker table or gene list)
# Outputs:  png (ComboPlot; FeatureDim and StackViolin are written next to it)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(qs2)
  library(tidyverse)
  library(Seurat)
  library(scplotter)
  library(plotthis)
  library(patchwork)
})
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "marker_features.R"))

p <- smk@params
out_base <- sub("__ComboPlot\\.png$", "", smk@output$png)
assay <- p$assay
res_col <- p$cluster_col
res_val <- sub(".*_res\\.", "", res_col)
title <- paste0("Resolution ", res_val, " (", res_col, ")")
subtitle <- toupper(gsub("_", " ", gsub("__", ": ", p$prefix)))

seu <- qs2::qs_read(smk@input$seurat_obj)
markers <- present_features(seu, resolve_features(p$spec, smk@input$features_file), assay)
log_message(sprintf("%d features: %s", length(markers), paste(markers, collapse = ", ")), "INFO")

if (length(markers) == 0) {
  placeholder_plot(paste0(out_base, "__ComboPlot"), p$prefix,
                   "No features to plot (no markers passed the filter, or none are in the assay)")
  log_message("No features; placeholder written", "DONE")
  quit(save = "no", status = 0)
}

n_markers <- length(markers)
few <- n_markers < 11

# Cluster labels padded with a space (as in the original), in cluster order
clusters <- seu@meta.data[[res_col]]
lvls <- if (is.factor(clusters)) levels(clusters) else as.character(sort(unique(clusters)))
seu@meta.data$Cluster <- factor(paste0(clusters, " "), levels = paste0(lvls, " "))
n_clust <- length(lvls)

# -- FEATURE UMAPS --------------------------------------------------------------
feature_dim <- function(ncol, palreverse) {
  FeatureStatPlot(
    seu,
    assay = assay,
    plot_type = "dim",
    features = markers,
    reduction = p$reduction,
    bg_cutoff = -Inf,
    palette = "ocean.dense",
    palreverse = palreverse,
    facet_ncol = ncol,
    raster = FALSE
  )
}

fd <- feature_dim(3, FALSE) +
  labs(title = subtitle) +
  theme(plot.title = element_text(size = 24), strip.text.x = element_text(size = 18))
plot_pdf_and_png(fd, paste0(out_base, "__FeatureDim"), width = 15,
                 height = max(4, 18 * n_markers / 10), png_only = TRUE)

# -- STACKED VIOLINS ------------------------------------------------------------
sv <- FeatureStatPlot(
  seu,
  assay = assay,
  title = title,
  subtitle = subtitle,
  plot_type = "violin",
  features = markers,
  ident = res_col,
  add_bg = TRUE,
  add_box = TRUE,
  box_ptsize = 1,
  stack = TRUE,
  palette = "jcolors-rainbow",
  palreverse = FALSE,
  raster = FALSE
) +
  theme(
    plot.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.subtitle = element_text(size = 12)
  ) +
  NoLegend()
plot_pdf_and_png(sv, paste0(out_base, "__StackViolin"), width = 3.5 + (n_clust / 2),
                 height = max(4, 14 * n_markers / 10), png_only = TRUE)

# -- COMBO PLOT -----------------------------------------------------------------
feat_ncol <- if (few) 5 else 12
row1_margin <- if (few) 60 else -10
row2_t_margin <- -30
row2_margin <- -40
pal <- "jcolors-rainbow"

p3 <- feature_dim(feat_ncol, TRUE) +
  guides(colour = guide_colorbar(title = "Expression")) +
  theme(
    plot.margin = margin(t = row1_margin, b = row1_margin),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, -3)
  )

p1 <- FeatureStatPlot(
  seu,
  assay = assay,
  plot_type = "violin",
  features = markers,
  ident = "Cluster",
  add_bg = TRUE,
  add_box = TRUE,
  box_ptsize = 1,
  stack = TRUE,
  flip = !few,
  palette = pal,
  palreverse = FALSE,
  raster = FALSE
) +
  theme(
    plot.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(t = row2_t_margin, b = row2_margin)
  ) +
  NoLegend()

p2 <- CellDimPlot(
  seu,
  group_by = "Cluster",
  reduction = p$reduction,
  label = TRUE,
  label_insitu = TRUE,
  label_size = 5,
  label_repel = TRUE,
  label_repulsion = 10,
  label_bg_r = 0.001,
  label_bg = "white",
  label_pt_size = 0.1,
  label_pt_color = "#00000000",
  label_segment_color = "#00000000",
  label_fg = "black",
  palette = pal,
  palreverse = FALSE,
  raster = FALSE
) +
  theme(
    plot.title = element_text(size = 14),
    plot.margin = margin(t = row2_t_margin, b = row2_t_margin, l = 20)
  )

h <- if (few) 14 else 3 + (n_markers / 2)
w <- if (few) 16 else 6 + n_markers
ratios <- if (few) c(6, 0.01, 5) else c(7 + n_markers / 2, -2, 6 + n_markers / 2)

combo <- (p1 | p2) / plot_spacer() / p3 +
  plot_layout(heights = ratios) +
  plot_annotation(
    title = title,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5)
    )
  )
plot_pdf_and_png(combo, paste0(out_base, "__ComboPlot"), width = w, height = h, png_only = TRUE)
log_message("Marker plots written", "DONE")
