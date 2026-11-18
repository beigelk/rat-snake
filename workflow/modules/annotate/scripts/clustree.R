###############################################################################
# Rule:     clustree
# Purpose:  How one population's clusters split as the resolution grows,
#           and a UMAP of the clusters at every resolution
#
# Description:
#   Wired from the exploratory visualization/clustree.R, which read a fixed
#   object and wrote to fixed folders. The clustree plot and the per-resolution
#   DimPlots are kept (same sizes); the marker tables and dot plots in that
#   script are now the findallmarkers, dotplot and markerplots rules.
#
# Inputs:   seurat_obj (compiled population object)
# Outputs:  png (clustree; a PDF is written next to it), umaps (folder of
#           DimPlots, one per resolution)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(clustree)
})

p <- smk@params
seu <- qs2::qs_read(smk@input$seurat_obj)
prefix <- paste0("cluster__", p$pop_name, "__", p$graph_name, "_res.")
res_cols <- grep(paste0("^", gsub(".", "\\.", prefix, fixed = TRUE)), colnames(seu@meta.data), value = TRUE)
log_message(sprintf("%d resolutions: %s", length(res_cols), paste(sub(prefix, "", res_cols, fixed = TRUE), collapse = ", ")), "INFO")

if (length(res_cols) > 1) {
  tree <- clustree(seu, prefix = prefix) +
    ggtitle(paste(p$project, p$pop_name)) +
    theme(legend.spacing = unit(0.1, "cm"), legend.key.size = unit(0.7, "lines"))
} else {
  tree <- ggplot() + theme_void() +
    labs(title = paste(p$project, p$pop_name), subtitle = "clustree needs two or more resolutions")
}
plot_pdf_and_png(tree, tools::file_path_sans_ext(smk@output$png), width = 12, height = 8)

dir.create(smk@output$umaps, recursive = TRUE, showWarnings = FALSE)
for (col in res_cols) {
  plot_pdf_and_png(
    DimPlot(seu, group.by = col, label = TRUE, reduction = p$reduction),
    file.path(smk@output$umaps, paste0(p$project, "__DimPlot__", col)),
    height = 5, width = 5
  )
}
log_message("clustree done", "DONE")
