
# -- SNAKEMAKE VARIABLES ------------------- #
# input
soupx_object    = snakemake@input[['soupx_object']]
doublet_id_file = snakemake@input[['doublet_id_file']]
adt_dsb_object  = snakemake@input[['adt_dsb_object']]

# params
r_library             = snakemake@params[['r_library']]
sample_name           = snakemake@params[['sample_name']]
rna_min_feat_cutoff   = snakemake@params[['rna_min_feat_cutoff']]
rna_max_feat_cutoff   = snakemake@params[['rna_max_feat_cutoff']]
rna_percent_mt_cutoff = snakemake@params[['rna_percent_mt_cutoff']]
adt_feat_min          = snakemake@params[['adt_feat_min']]
adt_isocount_max      = snakemake@params[['adt_isocount_max']]

# output
seurat_object        = snakemake@output[['seurat_object']]
rna_filt_cutoffs     = snakemake@output[['rna_filt_cutoffs']]
rna_sig_pcs          = snakemake@output[['rna_sig_pcs']]
rna_qcplots_prefilt  = snakemake@output[['rna_qcplots_prefilt']]
rna_umap_prefilt     = snakemake@output[['rna_umap_prefilt']]
rna_vln_prefilt      = snakemake@output[['rna_vln_prefilt']]
rna_scat_prefilt     = snakemake@output[['rna_scat_prefilt']]
rna_qcplots_postfilt = snakemake@output[['rna_qcplots_postfilt']]
rna_umap_postfilt    = snakemake@output[['rna_umap_postfilt']]
rna_vln_postfilt     = snakemake@output[['rna_vln_postfilt']]
rna_scat_postfilt    = snakemake@output[['rna_scat_postfilt']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)


# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# -- LOG INPUTS --------------------------- #
print(snakemake)

# -- LOAD PACKAGES ------------------- #
# analysis
library(Seurat)
library(dsb)

# utility
library(qs2)
library(tidyverse)
library(ggplot2)

# plotting
library(cowplot)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 40000 * 1024^2)

# -- SOURCE SCRIPTS ------------------------- #
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('utils.R', 'file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)


# -- FUNCTIONS ------------------------------ #
filter_rna_assay = function(
    seu.obj,
    min.feat.cutoff,
    max.feat.cutoff,
    percent.mt.cutoff
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
    file = paste0(rna_filt_cutoffs),
    quote = FALSE,
    row.names = FALSE,
    sep = ','
  )
  
  print(paste0('Filtering ', sample.name))
  
  # Filter the object
  seu.obj = subset(
    x = seu.obj,
    subset = nFeature_RNA > min.feat.cutoff &
             nFeature_RNA < max.feat.cutoff &
             percent_mt < percent.mt.cutoff
  )
  
  doublet.ids = read.table(doublet_id_file, header = TRUE)
  doublet.ids = doublet.ids[, 'doublet_ids']
  print(paste0('Number of doublets: ', length(doublet.ids)))
    if (length(doublet.ids) >= 1) {
      print(paste0('Removing ', length(doublet.ids), ' doublets.'))
      seu.obj <- seu.obj[, !colnames(seu.obj) %in% doublet.ids]
    }

  # Label the sample.name in the Seurat metadata
  sample.name_id = rep(sample.name, dim(seu.obj@assays$RNA$counts)[2])
  names(sample.name_id) = colnames(x = seu.obj)
  
  seu.obj = AddMetaData(
    object = seu.obj,
    metadata = sample.name_id,
    col.name = 'sample.nameID'
  )
  
  return(seu.obj)
  
}

get_sig_pcs = function(
    seu.obj
  ) {
  
  sample.name = seu.obj@project.name
  
  print('Calculating number of PCs.')
  
  # Determine percent of variation associated with each PC
  pct = seu.obj[['pca']]@stdev / sum(seu.obj[['pca']]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu = cumsum(pct)
  df = cbind(pct, cumu)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation
  # associated with the PC as less than 5
  co1 = which(cumu > 90 & pct < 5)[1]
  
  # Determine the difference between variation of PC and subsequent PC
  co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  pcs = min(co1, co2)
  
  print('Writing out final PC value.')
  write.csv(
    x = pcs,
    file = rna_sig_pcs,
    quote = FALSE,
    row.names = FALSE
  )
  
  return(pcs)
}

run_sctransform_etc = function(seu.obj) {
  
  if (is.null(seu.obj@meta.data$percent_mt)) {
    seu.obj[["percent_mt"]] = PercentageFeatureSet(seu.obj, pattern = "^MT-")
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
  seu.obj = RunPCA(seu.obj, npcs = 50)
  sig.pcs = get_sig_pcs(seu.obj = seu.obj)
  seu.obj = FindNeighbors(object = seu.obj, dims = 1:sig.pcs)
  seu.obj = FindClusters(object = seu.obj, resolution = 0.1)
  seu.obj = RunUMAP(
    object = seu.obj,
    dims = 1:sig.pcs,
    return.model = TRUE,
    n.neighbors = 15L,
    min.dist = 0.05
  )

  return(seu.obj)
  
}

plotting_qc = function(
    seu.obj,
    plot.name,
    plot.title,
    qc.plots,
    umap.plot,
    vln.plot,
    scat.plot
) {
  
  sample.name = seu.obj@project.name
  
  if (plot.name == 'rnafilt') {
    cutoff.tab = read.table(
      file = paste0(rna_filt_cutoffs),
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

get_isotype_binding = function(seurat.obj) {
  ## Create Isotype count matrix and filter out cells with high Isotype binding
  DefaultAssay(seurat.obj) = "ADT"
  
  adt_count_mtx = GetAssayData(object = seurat.obj, layer = "counts") %>%
    as.data.frame()
  
  iso_count_mtx = adt_count_mtx[grepl("Isotype", rownames(adt_count_mtx)), ]
  iso_count = colSums(iso_count_mtx)
  
  seurat.obj$IsoCount = iso_count
  
  return(seurat.obj)
  
}


filter_adt_assay = function(seurat.obj, min.feat.cutoff, max.iso.cutoff) {
  
  DefaultAssay(seurat.obj) = "ADT"
  
  # Filter cells
  seurat.obj = subset(
    seurat.obj,
    subset = nFeature_ADT > min.feat.cutoff & IsoCount < max.iso.cutoff
  )
  
  return(seurat.obj)
  
}

process_adt_assay = function(seurat.obj, norm.method = "CLR") {
  # ADT PRE-PROCESSING

  # CLR vs DSB normalization
  log_msg(sprintf("Normalization method selected: %s", norm.method), level =  'INFO')
  if (norm.method == "CLR") {
  log_msg(sprintf("Running NormalizeDate() in Seurat (norm.method=%s)", norm.method), level =  'INFO')
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
    log_msg(sprintf("NOT running normalization (norm.method=%s)", norm.method), level =  'INFO')

  }
    DefaultAssay(seurat.obj) = 'ADT'
    seurat.obj = FindVariableFeatures(seurat.obj, verbose = TRUE)

    prots = rownames(seurat.obj@assays$ADT)[!grepl('Isotype', rownames(seurat.obj@assays$ADT))]
    VariableFeatures(seurat.obj) = prots
    seurat.obj = ScaleData(seurat.obj, verbose = TRUE)
    seurat.obj = RunPCA(seurat.obj, features = prots, verbose = TRUE)
  
    names(seurat.obj@assays[["ADT"]]@layers) = c("counts", "data", "scale.data")
    
    return(seurat.obj)  
}


# -- MAIN ---------------------------------- #
# Load object 
s = qs2::qs_read(file = soupx_object)

# Make sure they are integers
# RNA params -----------------
rna_min_feat_cutoff   = as.integer(rna_min_feat_cutoff) 
rna_percent_mt_cutoff = as.integer(rna_percent_mt_cutoff)
rna_percent_mt_cutoff = as.integer(rna_percent_mt_cutoff)
rna_percent_mt_cutoff = as.integer(rna_percent_mt_cutoff)

# Get the median + 3*MAD if max_feat_cutoff == 'auto'
if (rna_max_feat_cutoff == 'auto') {
    rna_max_feat_cutoff = round(median(s$nFeature_RNA) + 3 * mad(s$nFeature_RNA))
} else {
    rna_max_feat_cutoff = as.integer(rna_max_feat_cutoff)
}

# ADT params -----------------
adt_feat_min     = as.integer(adt_feat_min)
adt_isocount_max = as.integer(adt_isocount_max)

# Normalize data with SCTransform
s = run_sctransform_etc(s)

plotting_qc(
  seu.obj    = s,
  plot.name  = 'prefilt',
  plot.title = 'Pre-filtered',
  qc.plots   = rna_qcplots_prefilt,
  umap.plot  = rna_umap_prefilt,
  vln.plot   = rna_vln_prefilt,
  scat.plot  = rna_scat_prefilt
)

# Filter RNA assay
s = filter_rna_assay(
  seu.obj           = s,
  min.feat.cutoff   = rna_min_feat_cutoff,
  max.feat.cutoff   = rna_max_feat_cutoff,
  percent.mt.cutoff = rna_percent_mt_cutoff
)

# Add DSB to seurat
s = add_dsb_seurat(s, adt_dsb_object)

s = get_isotype_binding(s)

s = filter_adt_assay(
  seurat.obj = s,
  min.feat.cutoff = adt_feat_min,
  max.iso.cutoff = adt_isocount_max
)

s = process_adt_assay(s)

s = run_sctransform_etc(s)

plotting_qc(
  seu.obj    = s,
  plot.name  = 'rnafilt',
  plot.title = 'Filtered',
  qc.plots   = rna_qcplots_postfilt,
  umap.plot  = rna_umap_postfilt,
  vln.plot   = rna_vln_postfilt,
  scat.plot  = rna_scat_postfilt
)

qs2::qs_save(
  object = s,
  file   = seurat_object
)
