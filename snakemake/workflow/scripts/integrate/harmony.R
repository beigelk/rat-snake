
# -- SNAKEMAKE VARIABLES ------------------- #
# input
merged_seurat = snakemake@input[['merged_seurat']]

# params
r_library         = snakemake@params[['r_library']]
adt_markers_file  = snakemake@params[['adt_markers_file']]
# wnn_umap_basename = snakemake@params[['wnn_umap_basename']]

# output
harmonized_seurat = snakemake@output[['harmonized_seurat']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)


# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# -- LOG INPUTS --------------------------- #
# print(snakemake)

# -- LOAD PACKAGES ------------------- #
# analysis
library(Seurat)
library(harmony)

# utility
library(qs2)
library(tidyverse)

# plot
library(ggplot2)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 200000 * 1024^2)

# source('Rscripts/utils.R', print.eval = FALSE)
print('functions')
process_multimodal_seurat = function(
    multimodal,
    adt.markers.file
) {

    # PREPPING DATA FOR MULTIMODAL OPERATIONS
    # -----------------------------------------------------------------------------------------------
    # Running PCA for the RNA assay
    # GEX
    DefaultAssay(multimodal) = 'RNA'
    multimodal = JoinLayers(multimodal, layers = c('counts', 'data', 'scale.data'))

    DefaultAssay(multimodal) = "RNA"
    multimodal = multimodal %>%
        NormalizeData(verbose = TRUE) %>%
        FindVariableFeatures(selection.method = 'vst', verbose = TRUE) %>%
        ScaleData(
            # vars.to.regress = c('percent_mt', 'nCount_RNA'),
            verbose = TRUE
        ) %>%
        RunPCA(
            assay = "RNA",
            npcs = 50,
            reduction.key = "rnaPC_", 
            reduction.name = "pca_rna"
        )

    # Running PCA for the ADT assay
    # ADT
    markers.oi = read.csv(adt.markers.file, header = FALSE)$V1

    DefaultAssay(multimodal) = 'ADT'
    multimodal = JoinLayers(multimodal, layers = c('counts', 'data', 'scale.data'))

    DefaultAssay(multimodal) = 'ADT'
    multimodal = multimodal %>%
        ScaleData() %>% 
        RunPCA(
            assay = "ADT",
            npcs = 15,
            reduction.key = "adtPC_", 
            reduction.name = "pca_adt",
            features = markers.oi
        )

    # -----------------------------------------------------------------------------------------------

    return(multimodal)
}

harmonize_assay = function(
    merged.seurat,
    assay.name
) {

    # Running Harmony
    # -----------------------------------------------------------------------------------------------
    # Manually set variable features of merged Seurat object
    DefaultAssay(merged.seurat) <- assay.name
    
    # Run Harmony to integrate samples
    # By default, harmony will look up the pca cell embeddings and use these to run harmony.
    # Therefore, it assumes that the Seurat object has these embeddings already precomputed.
    # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
    set.seed(42)
    harmonized.seurat = RunHarmony(
        merged.seurat,
        group.by.vars = c("orig.ident"), 
        reduction = paste0("pca_", tolower(assay.name)),
        reduction.save = paste0("harmony_", tolower(assay.name))
        # assay.use = "RNA",
        # plot.convergence = TRUE,
        # max_iter = 100,
        # kmeans_init_nstart = 20,
        # kmeans_init_iter_max = 1000
    )
  # -----------------------------------------------------------------------------------------------
  
  return(harmonized.seurat)
}

run_wnn = function(
    seurat.obj,
    rna.dims = ncol(Embeddings(seurat.obj, "harmony_rna")),
    adt.dims = ncol(Embeddings(seurat.obj, "harmony_adt")),
    k.nn.val = 20,
    prune.SNN.val = 1/15,
    wnn.name = NULL
) {

    # Run WNN on Harmony reduction
    # ----------------------------------------------------------------
    seurat.obj = FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_adt"),
        knn.graph.name = paste0("wknn_harmony", wnn.name),
        snn.graph.name = paste0("wsnn_harmony", wnn.name),
        weighted.nn.name = paste0("weighted.nn_harmony", wnn.name),
        dims.list = list(1:rna.dims, 1:adt.dims),
        modality.weight.name = c("RNA.weight", 'ADT.weight'),
        k.nn = k.nn.val,
        prune.SNN = prune.SNN.val,
        verbose = TRUE
    )
    # ----------------------------------------------------------------

    return(seurat.obj)
}

run_umap = function(
    seurat.obj,
    wnn.name = NULL
) {
    # Run UMAP for the weighted nearest neighbors (using Method 2)
    # ----------------------------------------------------------------
    seurat.obj = RunUMAP(
        seurat.obj,
        nn.name = paste0("weighted.nn_harmony", wnn.name),
        reduction.name = paste0("wnn.umap_harmony", wnn.name),
        reduction.key = "wnnUMAP_"
    )
    # ----------------------------------------------------------------

    return(seurat.obj)
}

###################################################################################

# -- Load Seurat object ---------------------------------------------------------- #
seu = qs2::qs_read(
    file = merged_seurat
)

# -- Scale, PCA based on ADT marker list ----------------------------------------- #
print('scale and pca')
seu = process_multimodal_seurat(seu, adt.markers.file = adt_markers_file)

# -- RNA: Run Harmony (via harmony package) -------------------------------------- #
print('harmony rna')
seu = harmonize_assay(
    merged.seurat = seu,
    assay.name = "RNA"
)

# -- ADT: Run Harmony (via harmony package) -------------------------------------- #
print('harmony adt')
seu = harmonize_assay(
    merged.seurat = seu,
    assay.name = "ADT"
)

# -- Find multimodal neighbors and calculate wnn --------------------------------- #
print('wnn default')
seu = run_wnn(
    seurat.obj = seu,
    wnn.name = "_default"
)

# -- Calculate UMAP from WNN ----------------------------------------------------- #
print('umap')
seu = run_umap(seu, wnn.name = "_default")

# -- Find multimodal neighbors and calculate wnn --------------------------------- #
print('wnn custom')
seu = run_wnn(
    seurat.obj = seu,
    rna.dims = 15,
    adt.dims = 15,
    k.nn.val = 10,
    prune.SNN.val = 1/30
)

# -- Save object ----------------------------------------------------------------- #
print('saving')
qs2::qs_save(
  object = seu,
  file = harmonized_seurat
)
print('saved')