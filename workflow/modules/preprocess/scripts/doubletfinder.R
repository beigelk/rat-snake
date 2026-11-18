###############################################################################
# Rule:     doubletfinder
# Purpose:  Call doublets in one sample with DoubletFinder
#
# Description:
#   Unchanged analysis: cells are pre-filtered (nFeature and percent
#   mitochondrial cutoffs; max nFeature "auto" = median + 3 MAD), the pK with
#   the highest BCmetric is chosen, the expected doublet rate follows 10x's
#   table, and DoubletFinder runs twice (with and without the homotypic
#   adjustment). A cell is a doublet if either run calls it one.
#   Fixes: the percent-mito cutoff was truncated to an integer (15.5 -> 15);
#   the mitochondrial pattern comes from preprocess.mt_pattern.
#   New: a one-row summary table for the report.
#
# Inputs:   snakemake@input$soupx_object
# Outputs:  doublet_id_file (barcodes), doublet_pkplot_pdf, doublet_qc_png,
#           summary (TSV: sample, cells_tested, doublets, doublet_rate_expected, pK)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "seurat_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(qs2)
  library(tidyverse)
  library(future)
  library(ggplot2)
  library(cowplot)
  library(rcartocolor)
})

p <- smk@params
plan(multicore, workers = as.integer(p$workers))


#' Expected multiplet rate from the number of cells (10x's table)
#'
#' @param seurat.obj Seurat object.
#'
#' @return Numeric rate, 0.004 for up to 500 cells up to 0.08 above 10,000.
get_multiplet_rate <- function(seurat.obj) {
  cell.count <- ncol(seurat.obj)
  breaks <- c(0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, Inf)
  rates  <- c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076, 0.08)
  rates[findInterval(cell.count, breaks, left.open = TRUE)]
}


#' Run DoubletFinder and add a final DF_Classification column
#'
#' @param seurat.obj SCTransformed Seurat object with PCA and clusters.
#' @param pN Numeric. DoubletFinder pN.
#' @param pk.pdf Character. PDF of BCmetric by pK.
#'
#' @return List: `seurat` (with DF_Classification), `pK`, `rate`.
run_doubletfinder <- function(seurat.obj, pN, pk.pdf) {
  sig_pcs <- get_sig_pcs(seurat.obj)
  doublet_rate <- get_multiplet_rate(seurat.obj)

  # pK identification (no ground truth)
  pdf(file = pk.pdf)
  sweep.res.list <- paramSweep(seurat.obj, PCs = 1:sig_pcs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  print(ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_point())
  dev.off()

  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk <- as.numeric(levels(bcmvn.max$pK))[bcmvn.max$pK]

  # Homotypic doublet proportion
  homotypic.prop <- modelHomotypic(seurat.obj@meta.data$seurat_clusters)
  nExp_poi <- round(doublet_rate * nrow(seurat.obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  for (n in c(nExp_poi, nExp_poi.adj)) {
    seurat.obj <- doubletFinder(seurat.obj, PCs = 1:sig_pcs, pN = pN, pK = optimal.pk,
                                nExp = n, reuse.pANN = NULL, sct = TRUE)
  }

  df.cols <- grep("DF.classification", colnames(seurat.obj@meta.data), value = TRUE)
  doublet.class <- seurat.obj@meta.data %>%
    select(all_of(df.cols)) %>%
    mutate(DF_Classification = case_when(
      if_any(all_of(df.cols), ~ .x == "Doublet") ~ "Doublet", .default = "Singlet"
    ))
  seurat.obj <- AddMetaData(seurat.obj, metadata = doublet.class$DF_Classification, col.name = "DF_Classification")
  list(seurat = seurat.obj, pK = optimal.pk, rate = doublet_rate)
}


#' UMAP and nFeature violin coloured by doublet call
#'
#' @param seurat.obj Seurat object with DF_Classification.
#' @param png.file Character. Output PNG.
plot_doublets <- function(seurat.obj, png.file) {
  dim.pal <- carto_pal(7, "ag_GrnYl")[c(1, 4)]
  # With fewer than 2 doublets no violin is drawn and the colours shift
  n_doublets <- sum(seurat.obj$DF_Classification == "Doublet")
  vln.pal <- if (n_doublets < 2) rev(dim.pal) else dim.pal
  png(file = png.file, width = 11, height = 8.5, res = 300, units = "in")
  umap <- DimPlot(seurat.obj, group.by = "DF_Classification", cols = dim.pal) + NoAxes()
  vln <- VlnPlot(seurat.obj, features = "nFeature_RNA", group.by = "DF_Classification",
                 cols = vln.pal, pt.size = 0.1)
  print(plot_grid(umap, vln, rel_widths = c(2.5, 1)))
  dev.off()
}


# -- MAIN ----------------------------------------------------------------------
seu <- qs2::qs_read(smk@input$soupx_object)

min_feat <- as.numeric(p$min_feat_cutoff)
pct_mt   <- as.numeric(p$percent_mt_cutoff)
max_feat <- if (identical(as.character(p$max_feat_cutoff), "auto")) {
  round(median(seu$nFeature_RNA) + 3 * mad(seu$nFeature_RNA))
} else {
  as.numeric(p$max_feat_cutoff)
}
log_message(sprintf("Pre-filter: %s < nFeature_RNA < %s, percent_mt < %s", min_feat, max_feat, pct_mt), "INFO")

seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = p$mt_pattern)
seu <- subset(seu, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent_mt < pct_mt)

res <- run_doubletfinder(seu, pN = as.numeric(p$pN), pk.pdf = smk@output$doublet_pkplot_pdf)
seu <- res$seurat
plot_doublets(seu, smk@output$doublet_qc_png)

doublet.ids <- colnames(seu)[seu$DF_Classification == "Doublet"]
write.table(doublet.ids, file = smk@output$doublet_id_file, quote = FALSE, sep = "",
            col.names = "doublet_ids", row.names = FALSE)
write.table(
  data.frame(sample = smk@wildcards$sample, cells_tested = ncol(seu), doublets = length(doublet.ids),
             doublet_rate_expected = res$rate, pK = res$pK),
  smk@output$summary, sep = "\t", quote = FALSE, row.names = FALSE
)
log_message(sprintf("%d doublets among %d cells", length(doublet.ids), ncol(seu)), "DONE")
