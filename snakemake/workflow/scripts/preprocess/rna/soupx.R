
# source("renv/activate.R")
# Activating environment modules: R/4.5.3
# # Bootstrapping renv 1.2.2 ---------------------------------------------------
# Warning: unable to access index for repository https://packagemanager.posit.co/cran/latest/src/contrib:
#   internet routines cannot be loaded
# - Downloading renv ... FAILED
# Warning message:
# failed to bootstrap renv: the project will not be loaded.
# Reason: failed to download:
# All download methods failed
# Use `renv::activate()` to re-initialize the project. 

# -- SNAKEMAKE VARIABLES ------------------- #
# input
filt_path = snakemake@input[['filt_path']]
raw_path  = snakemake@input[['raw_path']]

# params
r_library   = snakemake@params[['r_library']]
sample_name = snakemake@params[['sample_name']]
npcs        = snakemake@params[['npcs']]
resolution  = snakemake@params[['resolution']]

# output
soupx_object  = snakemake@output[['soupx_object']]
soupx_figures = snakemake@output[['soupx_figure']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)

# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")


Sys.getenv('RENV_PROJECT')
Sys.getenv('R_LIBS_USER')
Sys.getenv('R_ENVIRON_USER')
Sys.getenv('R_PROFILE_USER')

print(getwd())

print(.libPaths())
# -- LOAD PACKAGES -------------------------- #
# analysis
library(Seurat)
library(SoupX)

# utility
library(qs2)
library(stringr)

# plotting
library(ggplot2)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 1000 * 1024^2)


# -- FUNCTIONS ------------------------------ #
sct_to_clusters = function(seu, npcs, resolution) {
  seu = SCTransform(seu, method = "glmGamPoi")
  seu = RunPCA(object = seu, npcs = npcs)
  seu = RunUMAP(object = seu, dims = 1:npcs) 
  seu = FindNeighbors(object = seu, dims = 1:npcs)
  seu = FindClusters(object = seu, resolution = resolution)

  return(seu)
}

run_soupx = function(
  filt.path,
  raw.path,
  npcs.val,
  resolution.val,
  sample.name,
  soupx.object,
  soupx.figures
) {
  # Load filtered feature barcode matrix
  # This will have the Gene Expression and Antibody Capture
  filt = Read10X(data.dir = filt.path)
  raw = Read10X(data.dir = raw.path)
  
  # Create a Seurat object and run Seurat processing
  filt.seurat = CreateSeuratObject(counts = filt$`Gene Expression`)
  filt.seurat = sct_to_clusters(filt.seurat, npcs = npcs.val, resolution = resolution.val)

  pdf(file = paste0(soupx.figures))

  # Print UMAP before SoupX adjustment
  print(
    DimPlot(
      object = filt.seurat,
      reduction = "umap"
    ) +
    ggtitle(paste(sample.name, "PRE SoupX"))
  )
  
  # Identify names of clusters
  toc.clusters = filt.seurat@meta.data$seurat_clusters
  names(toc.clusters) = names(filt.seurat@active.ident)
  
  # SoupX call
  sc = SoupChannel(raw$`Gene Expression`, filt$`Gene Expression`)
  sc = setClusters(sc, toc.clusters)
  message("Starting autoEstCount")
  sc = autoEstCont(sc)
  soupx.out = adjustCounts(sc, roundToInt = TRUE)
  
  # Use the SoupX adjusted data
  sc = CreateSeuratObject(counts = soupx.out, project = sample.name)
  
  rownames(filt$`Antibody Capture`) = str_replace_all(rownames(filt$`Antibody Capture`), "_", "-")
  adt.assay = CreateAssay5Object(counts = filt$`Antibody Capture`)
  
  # Add this assay to the Seurat object
  sc[["ADT"]] = adt.assay
  
  # Re-process the data now that the SoupX adjustments have been applied
  sc = sct_to_clusters(sc, npcs = npcs.val, resolution = resolution.val)
  
  # Print UMAP after SoupX adjustment
  print(
    DimPlot(
      object = sc,
      reduction = "umap"
    ) +
    ggtitle(paste(sample.name, "POST SoupX")))
  dev.off()

  qs2::qs_save(
    object = sc,
    file = soupx.object
  )
}


# -- MAIN --------------------------- #
run_soupx(
  filt.path      = filt_path,
  raw.path       = raw_path,
  npcs.val       = as.integer(npcs),
  resolution.val = as.numeric(resolution),
  sample.name    = sample_name,
  soupx.object   = soupx_object,
  soupx.figures  = soupx_figures
)
