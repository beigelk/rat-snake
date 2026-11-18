###############################################################################
# Rule:     compile
# Purpose:  One object per population with every resolution's clusters and
#           the cell metadata
#
# Description:
#   Unchanged: old cluster columns are dropped, every resolution's cluster
#   table is merged by barcode and added, and cluster columns become factors
#   in numeric order.
#   Fix: the metadata table (sample metadata and clonotypes) was an input
#   but was never added; its columns are now added to the population's cells.
#   New: a table of cells per cluster at every resolution, for the report.
#
# Inputs:   target_obj (population object), cluster_results (one TSV per
#           resolution), metadata_obj (cell metadata table)
# Outputs:  compiled_obj (qs2), cluster_sizes (TSV: column, cluster, cells)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
})

seu <- qs2::qs_read(file = smk@input$target_obj)

cluster_files  <- unlist(smk@input$cluster_results)
cluster_tables <- lapply(cluster_files, read.delim, sep = "\t", colClasses = "character", check.names = FALSE)
cluster_cols   <- vapply(cluster_tables, function(t) colnames(t)[2], character(1))

# Drop old cluster columns, then add every resolution
for (col in cluster_cols) seu@meta.data[[col]] <- NULL
all_clusters <- Reduce(function(x, y) merge(x, y, by = "cell_barcode", all = TRUE), cluster_tables) %>%
  column_to_rownames("cell_barcode")
seu <- AddMetaData(seu, metadata = all_clusters[colnames(seu), , drop = FALSE])

# Cell metadata (sample metadata, clonotypes) not already in the object
meta <- qs2::qs_read(smk@input$metadata_obj)
new_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
if (length(new_cols) > 0) {
  seu <- AddMetaData(seu, metadata = meta[colnames(seu), new_cols, drop = FALSE])
  log_message(sprintf("Added metadata columns: %s", paste(new_cols, collapse = ", ")), "INFO")
}

# Cluster columns as factors in numeric order
for (col in cluster_cols) {
  vals <- as.character(seu@meta.data[[col]])
  seu@meta.data[[col]] <- factor(vals, levels = as.character(sort(as.numeric(unique(vals)))))
}

sizes <- map_dfr(cluster_cols, function(col) {
  seu@meta.data %>% count(cluster = .data[[col]], name = "cells") %>% mutate(column = col, .before = 1)
})
write_tsv(sizes, smk@output$cluster_sizes)
qs2::qs_save(seu, file = smk@output$compiled_obj)
log_message(sprintf("%s: %d cells, %d resolutions", smk@params$pop_name, ncol(seu), length(cluster_cols)), "DONE")
