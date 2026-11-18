###############################################################################
# File:     seurat_helpers.R
# Purpose:  Seurat steps used by more than one rat-snake script
#
# Description:
#   get_sig_pcs() was copied in doubletfinder.R and construct.R; run_wnn()
#   and run_wnn_umap() in harmony.R and subset.R (with different defaults,
#   which now come from the config: integrate.wnn and annotate.subset_wnn).
#
# Usage (in a rule's script, after snake_init()):
#   source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "seurat_helpers.R"))
###############################################################################

#' Number of significant PCs
#'
#' The smaller of: the first PC where cumulative variance passes 90% and
#' the PC's own share is under 5%; and the last PC after which the drop in
#' variance between consecutive PCs is over 0.1 percentage points.
#'
#' @param seu.obj Seurat object with a `pca` reduction.
#'
#' @return Integer number of PCs.
get_sig_pcs <- function(seu.obj) {
  pct  <- seu.obj[["pca"]]@stdev / sum(seu.obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1  <- which(cumu > 90 & pct < 5)[1]
  co2  <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  # The original returned NA (and the next step failed) when either rule
  # found nothing, e.g. with few PCs; now the other rule, or every PC, is used
  pcs  <- suppressWarnings(min(c(co1, co2), na.rm = TRUE))
  if (!is.finite(pcs)) {
    pcs <- length(pct)
    log_message(sprintf("No significant-PC cut-off found; using all %d PCs", pcs), "WARN")
  }
  log_message(sprintf("Significant PCs: %s", pcs), "INFO")
  pcs
}


#' Read a number written as a number or a fraction such as "1/15"
#'
#' @param x Numeric or character.
#'
#' @return Numeric.
as_fraction <- function(x) {
  if (is.numeric(x)) return(x)
  parts <- strsplit(trimws(as.character(x)), "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) return(as.numeric(parts[1]) / as.numeric(parts[2]))
  as.numeric(x)
}


#' Weighted nearest neighbours on the Harmony RNA and ADT embeddings
#'
#' Builds the graphs wknn_harmony / wsnn_harmony and the neighbour object
#' weighted.nn_harmony, with modality weights RNA.weight and ADT.weight.
#'
#' @param seurat.obj Seurat object with harmony_rna and harmony_adt reductions.
#' @param k.nn.val Integer. Neighbours (k.nn).
#' @param prune.SNN.val Numeric or "1/15"-style fraction. prune.SNN.
#' @param rna.dims,adt.dims Integer or NULL. Dimensions used; NULL = all.
#' @param snn.graph.name Character. SNN graph name (default wsnn_harmony).
#'
#' @return The Seurat object with the WNN graphs.
run_wnn <- function(seurat.obj, k.nn.val, prune.SNN.val, rna.dims = NULL, adt.dims = NULL,
                    snn.graph.name = "wsnn_harmony") {
  if (snn.graph.name != "wsnn_harmony") {
    stop(sprintf("Unknown graph name '%s'; only wsnn_harmony is built", snn.graph.name))
  }
  if (is.null(rna.dims)) rna.dims <- ncol(Embeddings(seurat.obj, "harmony_rna"))
  if (is.null(adt.dims)) adt.dims <- ncol(Embeddings(seurat.obj, "harmony_adt"))
  prune.SNN.val <- as_fraction(prune.SNN.val)
  log_message(sprintf("WNN: RNA dims 1:%s, ADT dims 1:%s, k.nn %s, prune.SNN %.4f",
                      rna.dims, adt.dims, k.nn.val, prune.SNN.val), "INFO")
  FindMultiModalNeighbors(
    object = seurat.obj,
    reduction.list = list("harmony_rna", "harmony_adt"),
    knn.graph.name = "wknn_harmony",
    snn.graph.name = snn.graph.name,
    weighted.nn.name = "weighted.nn_harmony",
    dims.list = list(1:rna.dims, 1:adt.dims),
    modality.weight.name = c("RNA.weight", "ADT.weight"),
    k.nn = k.nn.val,
    prune.SNN = prune.SNN.val,
    verbose = TRUE
  )
}


#' UMAP from the WNN neighbours (reduction wnn.umap_harmony)
#'
#' @param seurat.obj Seurat object after run_wnn().
#'
#' @return The Seurat object with the wnn.umap_harmony reduction.
run_wnn_umap <- function(seurat.obj) {
  RunUMAP(
    seurat.obj,
    nn.name = "weighted.nn_harmony",
    reduction.name = "wnn.umap_harmony",
    reduction.key = "wnnUMAP_"
  )
}
