###############################################################################
# Rule:     cluster
# Purpose:  Cluster one population at one resolution and save the clusters
#
# Description:
#   Unchanged: old clustering columns of the graph are removed, then
#   FindClusters (Louvain, algorithm 1) on the population's graph
#   (wsnn_harmony). Only the table of clusters is saved; the compile rule
#   gathers every resolution into the population's object.
#
# Inputs:   snakemake@input$seurat_obj   object the population is clustered on
# Params:   pop_name, graph_name, resolution (as written in file and column names)
# Outputs:  snakemake@output$results_tsv  cell_barcode + cluster__{pop}__{graph}_res.{res}
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
})

pop_name   <- smk@params$pop_name
graph_name <- smk@params$graph_name
res_str    <- smk@params$resolution
resolution <- as.numeric(res_str)

seu <- qs2::qs_read(file = smk@input$seurat_obj)
DefaultAssay(seu) <- "RNA"

# Remove earlier clustering columns of this graph
for (col in grep(paste0("^", graph_name, "_res"), colnames(seu@meta.data), value = TRUE)) {
  seu@meta.data[, col] <- NULL
}

log_message(sprintf("FindClusters on %s at resolution %s", graph_name, res_str), "INFO")
seu <- FindClusters(seu, algorithm = 1, graph.name = graph_name, resolution = resolution)

seurat_col  <- grep(paste0("^", graph_name, "_res"), colnames(seu@meta.data), value = TRUE)[1]
cluster_col <- paste0("cluster__", pop_name, "__", graph_name, "_res.", res_str)
cluster_results <- tibble(
  cell_barcode = colnames(seu),
  !!cluster_col := as.character(seu@meta.data[[seurat_col]])
)

readr::write_tsv(cluster_results, file = smk@output$results_tsv)
log_message(sprintf("%d clusters saved to %s", length(unique(cluster_results[[cluster_col]])),
                    smk@output$results_tsv), "DONE")
