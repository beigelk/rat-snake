###############################################################################
# Rule:     findallmarkers
# Purpose:  Markers of every cluster of one population, resolution and assay
#
# Description:
#   Unchanged: FindAllMarkers (positive markers only) on the cluster column;
#   cluster, gene, pct.1, pct.2, avg_log2FC, p_val and p_val_adj are kept.
#   Fix: when no markers are found the table still has its header (the
#   filter step and the plots read it), instead of an empty file.
#
# Inputs:   seurat_obj (compiled population object)
# Outputs:  raw_tsv (and an .xlsx copy next to it)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(openxlsx)
})

MARKER_COLS <- c("cluster", "gene", "pct.1", "pct.2", "avg_log2FC", "p_val", "p_val_adj")


#' Positive markers of every cluster
#'
#' @param seu_obj Seurat object.
#' @param meta_col Character. Cluster column.
#' @param assay_name Character. Assay to test (RNA or ADT).
#'
#' @return Data frame with MARKER_COLS (no rows when nothing is found or the test fails).
find_all_markers <- function(seu_obj, meta_col, assay_name) {
  markers <- tryCatch(
    FindAllMarkers(seu_obj, assay = assay_name, group.by = meta_col, only.pos = TRUE, verbose = TRUE),
    error = function(e) {
      log_message(paste("FindAllMarkers failed:", conditionMessage(e)), "WARN")
      NULL
    }
  )
  if (is.null(markers) || nrow(markers) == 0) {
    return(setNames(data.frame(matrix(ncol = length(MARKER_COLS), nrow = 0)), MARKER_COLS))
  }
  markers %>% remove_rownames() %>% dplyr::select(all_of(MARKER_COLS))
}


seu <- qs2::qs_read(smk@input$seurat_obj)
cluster_col <- smk@params$cluster_col
if (!cluster_col %in% colnames(seu@meta.data)) stop(sprintf("Column '%s' not found in metadata.", cluster_col))

Idents(seu) <- cluster_col
DefaultAssay(seu) <- smk@params$assay
markers <- find_all_markers(seu, meta_col = cluster_col, assay_name = smk@params$assay)

write.table(markers, smk@output$raw_tsv, row.names = FALSE, quote = FALSE, sep = "\t")
openxlsx::write.xlsx(markers, sub("\\.tsv$", ".xlsx", smk@output$raw_tsv), rowNames = FALSE)
log_message(sprintf("%d markers for %s (%s)", nrow(markers), cluster_col, smk@params$assay), "DONE")
