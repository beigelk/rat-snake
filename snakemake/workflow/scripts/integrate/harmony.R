
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
process_multimodal_seurat = function(multimodal, adt.markers.file) {

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

harmonize_rna_assay = function(merged.seurat) {

    # Running Harmony
    # -----------------------------------------------------------------------------------------------
    # Manually set variable features of merged Seurat object
    DefaultAssay(merged.seurat) <- 'RNA'
    # # VariableFeatures(merged.seurat) <- integ.features
    
    # # Calculate PCs using manually set variable features
    # merged.seurat <- RunPCA(merged.seurat, assay = "RNA", npcs = 50)
    
    # Run Harmony to integrate samples
    # By default, harmony will look up the pca cell embeddings and use these to run harmony.
    # Therefore, it assumes that the Seurat object has these embeddings already precomputed.
    # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
    set.seed(42)
    harmonized.seurat = RunHarmony(
        merged.seurat,
        group.by.vars = c("orig.ident"), 
        reduction = "pca_rna",
        reduction.save = "harmony_rna"
        # assay.use = "RNA",
        # plot.convergence = TRUE,
        # max_iter = 100,
        # kmeans_init_nstart = 20,
        # kmeans_init_iter_max = 1000
    )
  # -----------------------------------------------------------------------------------------------
  
  return(harmonized.seurat)
}

harmonize_adt_assay = function(merged.seurat) {
  
  # Manually set variable features of merged Seurat object
  DefaultAssay(merged.seurat) = 'ADT'
  
  # Running Harmony
  # -----------------------------------------------------------------------------------------------
  # Run Harmony to integrate samples
  # By default, harmony will look up the pca cell embeddings and use these to run harmony.
  # Therefore, it assumes that the Seurat object has these embeddings already precomputed.
  # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
    set.seed(42)
    harmonized.seurat = RunHarmony(
        merged.seurat,
        group.by.vars = c("orig.ident"), 
        reduction = "pca_adt",
        reduction.save = "harmony_adt"
        # assay.use = "ADT",
        # plot.convergence = TRUE,
        # max_iter = 100,
        # kmeans_init_nstart = 20,
        # kmeans_init_iter_max = 1000
    )
  
  return(harmonized.seurat)
}

run_wnn = function(seurat.obj) {
    rna_dims = ncol(Embeddings(seurat.obj, "harmony_rna"))
    adt_dims = ncol(Embeddings(seurat.obj, "harmony_adt"))
    # Run WNN on Harmony reduction
    seurat.obj = FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_adt"),
        knn.graph.name = "wknn_harmony",
        snn.graph.name = "wsnn_harmony",
        weighted.nn.name = "weighted.nn_harmony",
        dims.list = list(1:rna_dims, 1:adt_dims),
        modality.weight.name = c("RNA.weight", 'ADT.weight'),
        verbose = TRUE
    )
}

run_umap = function(seurat.obj) {
    # Run UMAP for the weighted nearest neighbors (using Method 2)
    # ----------------------------------------------------------------
    seurat.obj = RunUMAP(
        seurat.obj,
        nn.name = "weighted.nn_harmony",
        reduction.name = "wnn.umap_harmony",
        reduction.key = "wnnUMAP_"
    )
    # ----------------------------------------------------------------
    return(seurat.obj)
}

# plot_umaps = function(seurat.obj, anno.col1, anno.col2)
# {
#     # PLOT DIMPLOTS OF WNN
#     # ----------------------------------------------------------------
#     p1 = DimPlot(
#         object = seurat.obj,
#         group.by = anno.col1,
#         reduction = 'wnn.umap_harmony',
#         shuffle = TRUE) +
#     guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

#     plot_pdf_and_png(
#         plot = p1,
#         filename = wnn_umap_basename,
#         height = 10,
#         width = 10
#     )
#     # ----------------------------------------------------------------
# }

###################################################################################

seu = qs2::qs_read(
    file = merged_seurat
)

print('scale and pca')
# Scale, PCA
seu = process_multimodal_seurat(seu, adt.markers.file = adt_markers_file)

# Run Harmony (via harmony package)
print('harmony rna')
seu = harmonize_rna_assay(seu)
print('harmony adt')
seu = harmonize_adt_assay(seu)

# Find multimodal neighbors and calculate wnn
print('wnn')
seu = run_wnn(seu)

# Calculate UMAP from WNN
print('umap')
seu = run_umap(seu)
# plot_umaps(seu, anno.col1 = 'orig.ident', anno.col2 = 'seurat_clusters')
print('saving')
qs2::qs_save(
  object = seu,
  file = harmonized_seurat
)
print('saved')