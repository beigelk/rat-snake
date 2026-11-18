###############################################################################
# Rule:     soupx
# Purpose:  Remove ambient RNA with SoupX and build the sample's first Seurat
#           object (SoupX-corrected RNA + raw ADT counts)
#
# Description:
#   Unchanged analysis: SCTransform, PCA, UMAP and clusters on the filtered
#   matrix; SoupX (autoEstCont + adjustCounts) with those clusters; then the
#   corrected counts are re-processed the same way. ADT names have "_"
#   replaced by "-" as Seurat requires.
#   New: the estimated contamination fraction is written to a one-row table
#   for the report.
#
# Inputs:   snakemake@params$filtered_matrix, $raw_matrix (Cell Ranger multi)
# Outputs:  snakemake@output$soupx_object  qs2 Seurat object
#           snakemake@output$soupx_figure  PDF, UMAP before and after
#           snakemake@output$summary       TSV: sample, cells, rho
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(qs2)
  library(stringr)
  library(ggplot2)
})

sample_name <- smk@wildcards$sample
npcs        <- as.integer(smk@params$npcs)
resolution  <- as.numeric(smk@params$resolution)


#' SCTransform, PCA, UMAP, neighbours and clusters
#'
#' @param seu Seurat object.
#' @param npcs Integer. Number of PCs computed and used.
#' @param resolution Numeric. FindClusters resolution.
#'
#' @return The processed Seurat object.
sct_to_clusters <- function(seu, npcs, resolution) {
  seu <- SCTransform(seu, method = "glmGamPoi")
  seu <- RunPCA(object = seu, npcs = npcs)
  seu <- RunUMAP(object = seu, dims = 1:npcs)
  seu <- FindNeighbors(object = seu, dims = 1:npcs)
  seu <- FindClusters(object = seu, resolution = resolution)
  seu
}


#' Run SoupX on one sample and build its Seurat object
#'
#' @param filt.path,raw.path Character. Cell Ranger filtered and raw matrix folders.
#' @param npcs.val Integer. PCs for sct_to_clusters().
#' @param resolution.val Numeric. Clustering resolution for sct_to_clusters().
#' @param sample.name Character. Sample name (the object's project name).
#' @param soupx.figures Character. PDF of UMAPs before and after SoupX.
#'
#' @return List: `seurat` (corrected object with an ADT assay), `rho`
#'   (estimated contamination fraction) and `cells`.
run_soupx <- function(filt.path, raw.path, npcs.val, resolution.val, sample.name, soupx.figures) {
  # Filtered and raw matrices hold Gene Expression and Antibody Capture
  filt <- Read10X(data.dir = filt.path)
  raw  <- Read10X(data.dir = raw.path)

  filt.seurat <- CreateSeuratObject(counts = filt$`Gene Expression`)
  filt.seurat <- sct_to_clusters(filt.seurat, npcs = npcs.val, resolution = resolution.val)

  pdf(file = soupx.figures)
  print(DimPlot(object = filt.seurat, reduction = "umap") + ggtitle(paste(sample.name, "PRE SoupX")))

  toc.clusters <- filt.seurat@meta.data$seurat_clusters
  names(toc.clusters) <- names(filt.seurat@active.ident)

  sc <- SoupChannel(raw$`Gene Expression`, filt$`Gene Expression`)
  sc <- setClusters(sc, toc.clusters)
  log_message("Starting autoEstCont", "INFO")
  sc <- autoEstCont(sc)
  rho <- mean(sc$metaData$rho)
  soupx.out <- adjustCounts(sc, roundToInt = TRUE)

  seu <- CreateSeuratObject(counts = soupx.out, project = sample.name)
  rownames(filt$`Antibody Capture`) <- str_replace_all(rownames(filt$`Antibody Capture`), "_", "-")
  seu[["ADT"]] <- CreateAssay5Object(counts = filt$`Antibody Capture`)

  seu <- sct_to_clusters(seu, npcs = npcs.val, resolution = resolution.val)
  print(DimPlot(object = seu, reduction = "umap") + ggtitle(paste(sample.name, "POST SoupX")))
  dev.off()

  list(seurat = seu, rho = rho, cells = ncol(seu))
}


res <- run_soupx(
  filt.path      = smk@params$filtered_matrix,
  raw.path       = smk@params$raw_matrix,
  npcs.val       = npcs,
  resolution.val = resolution,
  sample.name    = sample_name,
  soupx.figures  = smk@output$soupx_figure
)

qs2::qs_save(object = res$seurat, file = smk@output$soupx_object)
write.table(
  data.frame(sample = sample_name, cells = res$cells, soupx_rho = round(res$rho, 4)),
  smk@output$summary, sep = "\t", quote = FALSE, row.names = FALSE
)
log_message(sprintf("SoupX done for %s: %d cells, contamination %.3f", sample_name, res$cells, res$rho), "DONE")
