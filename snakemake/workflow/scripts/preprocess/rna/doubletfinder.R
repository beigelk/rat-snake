
# -- SNAKEMAKE VARIABLES ------------------- #
# input
soupx_object = snakemake@input[['soupx_object']]

# params
r_library         = snakemake@params[['r_library']]
sample_name       = snakemake@params[['sample_name']]
min_feat_cutoff   = snakemake@params[['min_feat_cutoff']]
max_feat_cutoff   = snakemake@params[['max_feat_cutoff']]
percent_mt_cutoff = snakemake@params[['percent_mt_cutoff']]

# output
doublet_id_file    = snakemake@output[['doublet_id_file']]
doublet_pkplot_pdf = snakemake@output[['doublet_pkplot_pdf']]
doublet_qc_png     = snakemake@output[['doublet_qc_png']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)

# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# -- LOAD PACKAGES ------------------------- #
# analysis
library(Seurat)
library(DoubletFinder)

# utility
library(qs2)
library(tidyverse)
library(future)

# plotting
library(ggplot2)
library(cowplot)
library(rcartocolor)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 36000 * 1024^2)
plan(multicore, workers = 4)


# -- FUNCTIONS ------------------------------ #
get_multiplet_rate = function(seurat.obj)
{
  print('determing multiplet rate')
  cell.count = ncol(seurat.obj)
  
  doublet.predictor = NULL
  
  if (cell.count > 0 && cell.count <= 500)
  { doublet.predictor = 0.004 }
  if (cell.count > 500 && cell.count <= 1000)
  { doublet.predictor = 0.008 }
  if (cell.count > 1000 && cell.count <= 2000)
  { doublet.predictor = 0.016 }
  if (cell.count > 2000 && cell.count <= 3000)
  { doublet.predictor = 0.023 }
  if (cell.count > 3000 && cell.count <= 4000)
  { doublet.predictor = 0.031 }
  if (cell.count > 4000 && cell.count <= 5000)
  { doublet.predictor = 0.039 }
  if (cell.count > 5000 && cell.count <= 6000)
  { doublet.predictor = 0.046 }
  if (cell.count > 6000 && cell.count <= 7000)
  { doublet.predictor = 0.054 }
  if (cell.count > 7000 && cell.count <= 8000)
  { doublet.predictor = 0.061 }
  if (cell.count > 8000 && cell.count <= 9000)
  { doublet.predictor = 0.069 }
  if (cell.count > 9000 && cell.count <= 10000)
  { doublet.predictor = 0.076 }
  if (cell.count > 10000)
  { doublet.predictor = 0.08 }
  
  return(doublet.predictor)
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
  
  return(pcs)
}

run_doubletfinder = function(seurat.obj)
{  
  # Set the pN (default)
  pN = 0.25
  
  sample = seurat.obj@project.name
  
  # Read in the significant PCs
  sig_pcs = get_sig_pcs(seurat.obj)
  
  # Determine the multiplet rate based on 10X recommendations
  doublet_rate = get_multiplet_rate(seurat.obj)
  
  ## pK Identification (no ground-truth)
  # ----------------------------------------------------------------------------
  pdf(file = doublet_pkplot_pdf)

  sweep.res.list = paramSweep(seurat.obj, PCs = 1:sig_pcs, sct = TRUE)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn = find.pK(sweep.stats)

  # Plot a nicer pK plot with the actual axis labels
  print(ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_point())
  dev.off()
  # ----------------------------------------------------------------------------
  
  # Get pK at max BCmetric
  # ----------------------------------------------------------------------------
  bcmvn.max = bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk = as.numeric(levels(bcmvn.max$pK))[bcmvn.max$pK]
  # ----------------------------------------------------------------------------
  
  # Homotypic Doublet Proportion Estimate
  # ----------------------------------------------------------------------------
  annotations = seurat.obj@meta.data$seurat_clusters
  homotypic.prop = modelHomotypic(annotations)
  nExp_poi = round(doublet_rate * nrow(seurat.obj@meta.data))
  nExp_poi.adj = round(nExp_poi * (1 - homotypic.prop))
  # ----------------------------------------------------------------------------
  
  # Run DoubletFinder with varying classification stringencies
  # ----------------------------------------------------------------------------
  seurat.obj = doubletFinder(
    seurat.obj, PCs = 1:sig_pcs, pN = pN, pK = optimal.pk,
    nExp = nExp_poi,
    reuse.pANN = NULL,
    sct = TRUE
  )
  
  seurat.obj = doubletFinder(
    seurat.obj, PCs = 1:sig_pcs, pN = pN, pK = optimal.pk,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = TRUE
  )
  # ----------------------------------------------------------------------------
  
  write_doublet_ids(seurat.obj)
  
  return(seurat.obj)
  
}

write_doublet_ids = function(seurat.obj)
{  
  # Get the names of the doublet columns (usually 1, but could be >1)
  doublet.col.list = colnames(seurat.obj@meta.data)[grepl("DF.classification",
                                                           colnames(seurat.obj@meta.data))]
  print(doublet.col.list)
  # plot list for all plots at all DF.classifications
  plot_list = vector(mode = 'list')
  
  # If any of the DoubletFinder columns classify a cell as a doublet,
  # define it as a Doublet in the final Doublet_Classification
  doublet.class = seurat.obj@meta.data %>%
    select(all_of(doublet.col.list)) %>%
    mutate(
        DF_Classification = case_when(
            if_any(all_of(doublet.col.list), ~ .x == 'Doublet') ~ 'Doublet', .default = 'Singlet'
        )
    )
  
  # Add final Doublet classification to the Seurat object
  seurat.obj = AddMetaData(seurat.obj, metadata = doublet.class$DF_Classification, col.name = 'DF_Classification')

  # ------- Plotting -------
  # Colors for doublet plotting
  # colors from: https://github.com/Nowosad/rcartocolor
  dim.pal = carto_pal(7, "ag_GrnYl")[c(1, 4)]
  
  # Set the vln.pal differently based on if there are less than 2 doublets
  # When there's not enough cells to make a violin, the colors shift
  vln.pal = ''
  if (nrow(seurat.obj@meta.data %>% filter(DF_Classification == "Doublet")) < 2) {
    vln.pal = rev(dim.pal)
  } else {
    vln.pal = dim.pal
  }

  png(file = doublet_qc_png, width = 11, height = 8.5, res = 300, units = 'in')
  umap = DimPlot(seurat.obj, group.by = 'DF_Classification', cols = dim.pal) + NoAxes()
  vln = VlnPlot(seurat.obj, features = 'nFeature_RNA', group.by = 'DF_Classification', cols = vln.pal, pt.size = 0.1)
  print(plot_grid(umap, vln, rel_widths = c(2.5, 1)))
  dev.off()
  # ------------------------

  # Get cell barcodes from cells ID'd as doublets (from seurat object metadata)
  doublet.ids = rownames(seurat.obj@meta.data[seurat.obj@meta.data[, 'DF_Classification'] == "Doublet", ])

  # Single file of doublet IDs, write IDs to file
  print('Writing doublet IDs (cell barcodes) to file.')
  write.table(doublet.ids, file = doublet_id_file, quote = FALSE, sep = '', col.names = 'doublet_ids', row.names = FALSE)
  print('Doublet IDs written.')
}


# -- MAIN ------------------------------------------------- #
seu = qs2::qs_read(soupx_object)

# Make sure they are integers
min_feat_cutoff   = as.integer(min_feat_cutoff)
percent_mt_cutoff = as.integer(percent_mt_cutoff)

# Get the median + 3*MAD if max_feat_cutoff == 'auto'
if (max_feat_cutoff == 'auto') {
    max_feat_cutoff = round(median(seu$nFeature_RNA) + 3 * mad(seu$nFeature_RNA))
} else {
    max_feat_cutoff = as.integer(max_feat_cutoff)
}

# Load Seurat object

# Filter
seu[["percent_mt"]] = PercentageFeatureSet(seu, pattern = "^MT-")
seu = subset(
    seu,
    subset = nFeature_RNA > min_feat_cutoff &
              nFeature_RNA < max_feat_cutoff &
              percent_mt < percent_mt_cutoff
)

# Run DoubletFinder
seu = run_doubletfinder(seu)
