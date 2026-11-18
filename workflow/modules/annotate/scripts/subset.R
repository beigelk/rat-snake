###############################################################################
# Rule:     subset
# Purpose:  Cut a population out of its parent's clusters and re-run WNN + UMAP
#
# Description:
#   Unchanged: the parent's clusters at base_resolution are added, cells in
#   base_clusters are kept ("all" keeps every cell), and FindMultiModalNeighbors
#   and the WNN UMAP are re-run on the Harmony embeddings.
#   Fixes: the cluster table is matched to cells by barcode (AddMetaData was
#   given a table without row names) and a mismatch now says which cells
#   differ (if (all.equal(...)) errored instead); the "cells after
#   subsetting" log line printed no number. WNN settings come from
#   annotate.subset_wnn (k_nn 10, prune_snn 1/30, 15 RNA and 15 ADT dims, as
#   before) or the population's own wnn block.
#
# Inputs:   seurat_obj (parent's object), cluster_results (parent's clusters)
# Outputs:  seurat_sub_obj (qs2)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "seurat_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
})

p <- smk@params
spec <- p$subset

seu <- qs2::qs_read(file = smk@input$seurat_obj)
DefaultAssay(seu) <- "RNA"

clusters <- readr::read_tsv(smk@input$cluster_results, col_types = cols(.default = col_character()))
missing <- setdiff(colnames(seu), clusters$cell_barcode)
if (length(missing) > 0) {
  stop(sprintf("%d of %d cells have no cluster in %s (e.g. %s)", length(missing), ncol(seu),
               smk@input$cluster_results, paste(head(missing, 3), collapse = ", ")))
}
seu <- AddMetaData(seu, metadata = column_to_rownames(clusters, "cell_barcode")[colnames(seu), , drop = FALSE])

cluster_col <- paste0("cluster__", spec$base_name, "__", spec$base_graph_name, "_res.", spec$base_resolution)
keep <- unlist(spec$base_clusters)
log_message(sprintf("Cells before subsetting: %d", ncol(seu)), "INFO")
if (identical(keep, "all")) {
  log_message(sprintf("base_clusters is 'all': keeping every cluster of %s", cluster_col), "INFO")
} else {
  log_message(sprintf("Keeping clusters %s of %s", paste(keep, collapse = ", "), cluster_col), "INFO")
  seu <- seu[, seu@meta.data[[cluster_col]] %in% as.character(keep)]
}
log_message(sprintf("Cells after subsetting: %d", ncol(seu)), "INFO")

wnn <- p$wnn
seu <- run_wnn(seu, k.nn.val = wnn$k_nn, prune.SNN.val = wnn$prune_snn,
               rna.dims = wnn$rna_dims, adt.dims = wnn$adt_dims, snn.graph.name = p$graph_name)
seu <- run_wnn_umap(seu)

qs2::qs_save(seu, file = smk@output$seurat_sub_obj)
log_message("Subset saved", "DONE")
