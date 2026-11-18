###############################################################################
# Rule:     construct
# Purpose:  Build one sample's filtered multimodal Seurat object (RNA + ADT)
#
# Description:
#   Unchanged analysis: SCTransform + PCA/UMAP of the SoupX object (QC plots
#   before filtering); RNA filter (nFeature bounds, max nFeature "auto" =
#   median + 3 MAD, percent mito) and doublet removal; dsb values into the
#   ADT data layer; IsoCount and ADT filter; ADT normalisation (CLR by
#   default, as before), scaling and PCA; SCTransform again and QC plots
#   after filtering.
#   Fixes: the percent-mito cutoff is no longer truncated to an integer; the
#   sample column is Sample_ID (was sample.nameID); a Dataset column is added.
#   New: preprocess.construct.adt.normalization (CLR or dsb), npcs and
#   cluster_resolution options, and a summary table of cells kept per step.
#
# Inputs:   soupx_object, doublet_id_file, adt_dsb_object
# Outputs:  seurat_object (qs2), rna_filt_cutoffs, rna_sig_pcs, QC plots
#           before (prefilt) and after (postfilt), summary (TSV)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "seurat_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
})

p   <- smk@params
out <- smk@output

# -- FUNCTIONS -----------------------------------------------------------------
#' Filter cells on RNA QC and remove doublets
#'
#' Writes the cutoffs used to the rna_filt_cutoffs table and adds the
#' Sample_ID column (the sample name).
#'
#' @param seu.obj Seurat object.
#' @param min.feat.cutoff,max.feat.cutoff Numeric. nFeature_RNA bounds (exclusive).
#' @param percent.mt.cutoff Numeric. Maximum percent mitochondrial counts.
#' @param doublet.file Character. File of doublet barcodes from doubletfinder.
#' @param cutoffs.file Character. Output table of the cutoffs.
#'
#' @return The filtered Seurat object.
filter_rna_assay = function(
    seu.obj,
    min.feat.cutoff,
    max.feat.cutoff,
    percent.mt.cutoff,
    doublet.file,
    cutoffs.file
  ) {
  
  sample.name = seu.obj@project.name
  DefaultAssay(seu.obj) = "RNA"
  
  # Filter cells
  cutoffs = data.frame(
    nFeature_RNA_min = min.feat.cutoff,
    nFeature_RNA_max = max.feat.cutoff,
    percent_mt_max = percent.mt.cutoff
  )
  
  write.table(
    cutoffs,
    file = cutoffs.file,
    quote = FALSE,
    row.names = FALSE,
    sep = ','
  )
  
  log_message(paste0('Filtering ', sample.name), 'INFO')
  
  # Filter the object
  seu.obj = subset(
    x = seu.obj,
    subset = nFeature_RNA > min.feat.cutoff &
             nFeature_RNA < max.feat.cutoff &
             percent_mt < percent.mt.cutoff
  )
  
  doublet.ids = read.table(doublet.file, header = TRUE)
  doublet.ids = doublet.ids[, 'doublet_ids']
  log_message(paste0('Number of doublets: ', length(doublet.ids)), 'INFO')
    if (length(doublet.ids) >= 1) {
      log_message(paste0('Removing ', length(doublet.ids), ' doublets.'), 'INFO')
      seu.obj <- seu.obj[, !colnames(seu.obj) %in% doublet.ids]
    }

  # Label the sample in the Seurat metadata (was the column sample.nameID)
  seu.obj$Sample_ID = sample.name
  
  return(seu.obj)
  
}

#' SCTransform (regressing percent_mt), PCA, neighbours, clusters and UMAP
#'
#' @param seu.obj Seurat object.
#' @param npcs Integer. PCs computed; the significant ones (get_sig_pcs) are used.
#' @param resolution Numeric. FindClusters resolution.
#' @param mt.pattern Character. Regex for mitochondrial genes.
#' @param sig.pcs.file Character. File the number of significant PCs is written to.
#'
#' @return The processed Seurat object.
run_sctransform_etc = function(seu.obj, npcs, resolution, mt.pattern, sig.pcs.file) {
  
  if (is.null(seu.obj@meta.data$percent_mt)) {
    seu.obj[["percent_mt"]] = PercentageFeatureSet(seu.obj, pattern = mt.pattern)
  }

  ## Normalization
  # all.genes = rownames(seu.obj)
  seu.obj = SCTransform(
    object = seu.obj,
    vars.to.regress = c('percent_mt'),
    return.only.var.genes = FALSE,
    verbose = TRUE
  )

  ## Perform linear dimensional reduction
  seu.obj = RunPCA(seu.obj, npcs = npcs)
  sig.pcs = get_sig_pcs(seu.obj = seu.obj)
  write.csv(x = sig.pcs, file = sig.pcs.file, quote = FALSE, row.names = FALSE)
  seu.obj = FindNeighbors(object = seu.obj, dims = 1:sig.pcs)
  seu.obj = FindClusters(object = seu.obj, resolution = resolution)
  seu.obj = RunUMAP(
    object = seu.obj,
    dims = 1:sig.pcs,
    return.model = TRUE,
    n.neighbors = 15L,
    min.dist = 0.05
  )

  return(seu.obj)
  
}

#' QC plots: UMAP by cluster, violins and feature scatters
#'
#' @param seu.obj Seurat object.
#' @param plot.name Character. "rnafilt" adds the cutoffs as a subtitle.
#' @param plot.title Character. Title after the sample name.
#' @param qc.plots Character. PDF with all three plots.
#' @param umap.plot,vln.plot,scat.plot Character. PNG of each plot.
#' @param cutoffs.file Character. Table written by filter_rna_assay().
plotting_qc = function(
    seu.obj,
    plot.name,
    plot.title,
    qc.plots,
    umap.plot,
    vln.plot,
    scat.plot,
    cutoffs.file = NULL
) {
  
  sample.name = seu.obj@project.name
  
  if (plot.name == 'rnafilt') {
    cutoff.tab = read.table(
      file = cutoffs.file,
      header = TRUE, sep = ','
    )

    cutoff_info = paste0(
      'nFeature_RNA min: ', cutoff.tab$nFeature_RNA_min, '; ',
      'nFeature_RNA max (3xMAD): ', cutoff.tab$nFeature_RNA_max, '; ',
      'percent_mito max: ', cutoff.tab$percent_mt_max
    )
  } else {
    cutoff_info = 'No filtering'
  }
  
  num_cells = ncol(seu.obj)
  
  # Making plots
  umap_clust = DimPlot(
    seu.obj,
    group.by = "seurat_clusters"
  ) +
  labs(
    title = paste(sample.name, plot.title, sep = ' '),
    subtitle = cutoff_info, 
    caption = paste('Number of cells:', num_cells))
  
  vln_rna = plot_grid(
    ggplot() +
      labs(
        title = paste(sample.name, plot.title, sep = ' '),
        subtitle = cutoff_info, 
        caption = paste('Number of cells:', num_cells)) +
      theme(plot.title = element_text(size = 15, vjust = -0.5, face = 'bold')),
    VlnPlot(seu.obj, features = c("nCount_RNA", "nFeature_RNA", "percent_mt")),
    ncol = 1,
    rel_heights = c(0.1, 1))
  
  feat_scat1 = FeatureScatter(seu.obj,
                              feature1 = 'nCount_RNA', feature2 = 'percent_mt',
                              pt.size = 0.3)
  
  feat_scat2 = FeatureScatter(seu.obj,
                              feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA',
                              pt.size = 0.3)
  
  feat_scat = plot_grid(
    ggplot() +
      labs(title =
       paste(sample.name,
       plot.title, sep = ' '),
           subtitle = cutoff_info, 
           caption = paste('Number of cells:', num_cells)) +
      theme(plot.title = element_text(size = 15, vjust = -0.5, face = 'bold')),
    (feat_scat1 | feat_scat2),
    ncol = 1,
    rel_heights = c(0.1, 1))
  
  # Plotting to file
  pdf(file = paste0(qc.plots),
      width = 8, height = 8)
  print(umap_clust)
  print(vln_rna)
  print(feat_scat)
  dev.off()
  
  png(file = paste0(umap.plot),
      width = 8, height = 8, res = 300, units = "in")
  print(umap_clust)
  dev.off()

  png(file = paste0(vln.plot),
      width = 12, height = 8, res = 300, units = "in")
  print(vln_rna)
  dev.off()
  
  png(file = paste0(scat.plot),
      width = 8, height = 8, res = 300, units = "in")
  print(feat_scat)
  dev.off()
  
}

#' Put the dsb-normalised ADT values in the ADT assay's data layer
#'
#' Keeps only cells present in both the object and the dsb output.
#'
#' @param seurat.obj Seurat object with an ADT assay.
#' @param adt.dsb Character. qs2 file from the dsb rule.
#'
#' @return The Seurat object restricted to cells dsb kept.
add_dsb_seurat = function(seurat.obj, adt.dsb) {
    DefaultAssay(seurat.obj) = 'ADT'
    sample = seurat.obj@project.name
    dsb.out = qs2::qs_read(file = adt.dsb)

    # Get the common cells between the Seurat object and the dsb data
    cells = intersect(
      colnames(seurat.obj),
      colnames(dsb.out$cells.dsb.norm$dsb_normalized_matrix)
    )

    dsb.norm.mat = dsb.out$cells.dsb.norm$dsb_normalized_matrix[, cells]
    cells.to.keep = intersect(colnames(seurat.obj), colnames(dsb.norm.mat))

    # print(cells.to.keep)

    seurat.obj_subset = subset(seurat.obj, cells = cells.to.keep)
    dsb.norm.mat.new = dsb.norm.mat[, cells.to.keep]

    rownames(dsb.norm.mat.new) = gsub('_', '-', rownames(dsb.norm.mat.new))

    seurat.obj.final = SetAssayData(
        seurat.obj_subset,
        layer = "data",
        dsb.norm.mat.new
    )

    return(seurat.obj.final)
  
}

#' Add IsoCount: total isotype control counts per cell
#'
#' @param seurat.obj Seurat object with an ADT assay.
#' @param isotype.pattern Character. Regex for isotype control names.
#'
#' @return The Seurat object with an IsoCount column.
get_isotype_binding = function(seurat.obj, isotype.pattern = "Isotype") {
  ## Create Isotype count matrix and filter out cells with high Isotype binding
  DefaultAssay(seurat.obj) = "ADT"
  
  adt_count_mtx = GetAssayData(object = seurat.obj, layer = "counts") %>%
    as.data.frame()
  
  iso_count_mtx = adt_count_mtx[grepl(isotype.pattern, rownames(adt_count_mtx)), ]
  iso_count = colSums(iso_count_mtx)
  
  seurat.obj$IsoCount = iso_count
  
  return(seurat.obj)
  
}


#' Filter cells on ADT QC
#'
#' @param seurat.obj Seurat object with nFeature_ADT and IsoCount.
#' @param min.feat.cutoff Numeric. Minimum nFeature_ADT (exclusive).
#' @param max.iso.cutoff Numeric. Maximum IsoCount (exclusive).
#'
#' @return The filtered Seurat object.
filter_adt_assay = function(seurat.obj, min.feat.cutoff, max.iso.cutoff) {
  
  DefaultAssay(seurat.obj) = "ADT"
  
  # Filter cells
  seurat.obj = subset(
    seurat.obj,
    subset = nFeature_ADT > min.feat.cutoff & IsoCount < max.iso.cutoff
  )
  
  return(seurat.obj)
  
}

#' Normalise (CLR) or keep dsb values, then scale and PCA the ADT assay
#'
#' Isotype controls are left out of the variable features.
#'
#' @param seurat.obj Seurat object.
#' @param norm.method Character. "CLR" (original behaviour; replaces the dsb
#'   values in the data layer) or "dsb" (keeps them).
#' @param isotype.pattern Character. Regex for isotype control names.
#'
#' @return The Seurat object with a PCA of the ADT assay.
process_adt_assay = function(seurat.obj, norm.method = "CLR", isotype.pattern = "Isotype") {
  # ADT PRE-PROCESSING

  # CLR vs DSB normalization
  log_message(sprintf("Normalization method selected: %s", norm.method), 'INFO')
  if (norm.method == "CLR") {
  log_message(sprintf("Running NormalizeData() in Seurat (norm.method=%s)", norm.method), 'INFO')
    seurat.obj <- NormalizeData(
      seurat.obj,
      assay = "ADT",
      normalization.method = "CLR",
      margin = 2
    )
  } else if (norm.method == "dsb") {
    # Seurat multimodal with adaptations for DSB-normalized data
    # notice the normalization is excluded here, that's because the data has
    # already been normalized with DSB.
    # https://github.com/satijalab/seurat/issues/3890
    # !!! Do not normalize here (dsb did that) !!!
    log_message(sprintf("NOT running normalization (norm.method=%s)", norm.method), 'INFO')

  }
    DefaultAssay(seurat.obj) = 'ADT'
    seurat.obj = FindVariableFeatures(seurat.obj, verbose = TRUE)

    prots = rownames(seurat.obj@assays$ADT)[!grepl(isotype.pattern, rownames(seurat.obj@assays$ADT))]
    VariableFeatures(seurat.obj) = prots
    seurat.obj = ScaleData(seurat.obj, verbose = TRUE)
    seurat.obj = RunPCA(seurat.obj, features = prots, verbose = TRUE)
  
    names(seurat.obj@assays[["ADT"]]@layers) = c("counts", "data", "scale.data")
    
    return(seurat.obj)  
}


# -- MAIN ----------------------------------------------------------------------
s <- qs2::qs_read(file = smk@input$soupx_object)
n_start <- ncol(s)

rna_min_feat_cutoff   <- as.numeric(p$rna$rna_min_feat_cutoff)
rna_percent_mt_cutoff <- as.numeric(p$rna$rna_percent_mt_cutoff)
rna_max_feat_cutoff   <- if (identical(as.character(p$rna$rna_max_feat_cutoff), "auto")) {
  round(median(s$nFeature_RNA) + 3 * mad(s$nFeature_RNA))
} else {
  as.numeric(p$rna$rna_max_feat_cutoff)
}
adt_feat_min     <- as.numeric(p$adt$adt_feat_min)
adt_isocount_max <- as.numeric(p$adt$adt_isocount_max)
iso_pattern      <- p$isotype_pattern

sct_args <- list(npcs = as.integer(p$npcs), resolution = as.numeric(p$cluster_resolution),
                 mt.pattern = p$mt_pattern, sig.pcs.file = out$rna_sig_pcs)

s <- do.call(run_sctransform_etc, c(list(seu.obj = s), sct_args))
plotting_qc(s, plot.name = "prefilt", plot.title = "Pre-filtered",
            qc.plots = out$rna_qcplots_prefilt, umap.plot = out$rna_umap_prefilt,
            vln.plot = out$rna_vln_prefilt, scat.plot = out$rna_scat_prefilt)

s <- filter_rna_assay(s, min.feat.cutoff = rna_min_feat_cutoff, max.feat.cutoff = rna_max_feat_cutoff,
                      percent.mt.cutoff = rna_percent_mt_cutoff, doublet.file = smk@input$doublet_id_file,
                      cutoffs.file = out$rna_filt_cutoffs)
n_rna <- ncol(s)

s <- add_dsb_seurat(s, smk@input$adt_dsb_object)
n_dsb <- ncol(s)

s <- get_isotype_binding(s, isotype.pattern = iso_pattern)
s <- filter_adt_assay(s, min.feat.cutoff = adt_feat_min, max.iso.cutoff = adt_isocount_max)
n_adt <- ncol(s)

s <- process_adt_assay(s, norm.method = p$adt$normalization, isotype.pattern = iso_pattern)
s <- do.call(run_sctransform_etc, c(list(seu.obj = s), sct_args))
s$Dataset <- smk@wildcards$dataset

plotting_qc(s, plot.name = "rnafilt", plot.title = "Filtered",
            qc.plots = out$rna_qcplots_postfilt, umap.plot = out$rna_umap_postfilt,
            vln.plot = out$rna_vln_postfilt, scat.plot = out$rna_scat_postfilt,
            cutoffs.file = out$rna_filt_cutoffs)

qs2::qs_save(object = s, file = out$seurat_object)
write.table(
  data.frame(sample = smk@wildcards$sample, dataset = smk@wildcards$dataset,
             cells_soupx = n_start, cells_after_rna_filter = n_rna, cells_after_dsb = n_dsb,
             cells_final = n_adt, nFeature_RNA_min = rna_min_feat_cutoff, nFeature_RNA_max = rna_max_feat_cutoff,
             percent_mt_max = rna_percent_mt_cutoff, adt_feat_min = adt_feat_min, adt_isocount_max = adt_isocount_max),
  out$summary, sep = "\t", quote = FALSE, row.names = FALSE
)
log_message(sprintf("%s: %d -> %d cells", smk@wildcards$sample, n_start, n_adt), "DONE")
