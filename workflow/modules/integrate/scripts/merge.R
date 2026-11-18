###############################################################################
# Rule:     merge
# Purpose:  Merge every sample's Seurat object into one
#
# Description:
#   Unchanged: objects are named by their project name (the sample) and
#   merged with merge.data = TRUE; cell names get the sample as a prefix
#   ("{sample}_{barcode}"), which is also how scRepertoire names cells.
#
# Inputs:   snakemake@input$seurat_objects  one qs2 per sample
# Outputs:  snakemake@output$merged_seurat  qs2
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
})


#' Read Seurat objects into a list named by each object's project name
#'
#' @param seurat.file.list Character vector of qs2 files.
#'
#' @return Named list of Seurat objects.
make_seurat_list <- function(seurat.file.list) {
  seurat.obj.list <- list()
  for (filename in seurat.file.list) {
    seu <- qs2::qs_read(file = filename)
    log_message(sprintf("%s: %d cells", seu@project.name, ncol(seu)), "INFO")
    seurat.obj.list[[seu@project.name]] <- seu
  }
  seurat.obj.list
}


#' Merge a list of Seurat objects (Seurat 5), prefixing cell names with the list names
#'
#' @param seurat.obj.list Named list of Seurat objects.
#'
#' @return The merged Seurat object (the only object when the list has one).
make_merged_seurat <- function(seurat.obj.list) {
  if (length(seurat.obj.list) == 1) {
    return(RenameCells(seurat.obj.list[[1]], add.cell.id = names(seurat.obj.list)[1]))
  }
  merge(
    x = seurat.obj.list[[1]],
    y = seurat.obj.list[2:length(seurat.obj.list)],
    merge.data = TRUE,
    add.cell.ids = names(seurat.obj.list)
  )
}


merged <- make_merged_seurat(make_seurat_list(unlist(smk@input$seurat_objects)))
qs2::qs_save(object = merged, file = smk@output$merged_seurat)
log_message(sprintf("Merged %d cells", ncol(merged)), "DONE")
