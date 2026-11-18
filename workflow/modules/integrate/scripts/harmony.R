###############################################################################
# Rule:     harmony
# Purpose:  Integrate the samples: Harmony on the RNA and ADT PCAs, then
#           weighted nearest neighbours (WNN) and a WNN UMAP
#
# Description:
#   Unchanged analysis. Settings that were written into the script are now
#   options with the same defaults: integrate.harmony (group_by_vars
#   [orig.ident], rna_npcs 50, adt_npcs 15, adt_markers_file) and
#   integrate.wnn (k_nn 20, prune_snn 1/15, all Harmony dimensions).
#   Without adt_markers_file the ADT PCA uses every non-isotype ADT.
#
# Inputs:   snakemake@input$merged_seurat
# Outputs:  snakemake@output$harmonized_seurat (graphs wknn_harmony and
#           wsnn_harmony, reduction wnn.umap_harmony)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)
source(file.path(smk@scriptdir, "..", "..", "..", "lib", "R", "seurat_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(qs2)
  library(tidyverse)
})

h    <- smk@params$harmony
wnn  <- smk@params$wnn
seed <- if (is.null(smk@config$seed)) 42 else smk@config$seed

# -- FUNCTIONS -----------------------------------------------------------------
#' Normalise, scale and PCA the RNA and ADT assays of the merged object
#'
#' RNA: NormalizeData, FindVariableFeatures (vst), ScaleData, PCA (pca_rna).
#' ADT: ScaleData and PCA (pca_adt) on the listed markers.
#'
#' @param multimodal Merged Seurat object.
#' @param adt.markers Character vector of ADT features for the PCA.
#' @param rna.npcs,adt.npcs Integer. PCs per assay (50 and 15 before).
#'
#' @return The Seurat object with pca_rna and pca_adt.
process_multimodal_seurat = function(multimodal, adt.markers, rna.npcs, adt.npcs) {

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
            npcs = rna.npcs,
            reduction.key = "rnaPC_", 
            reduction.name = "pca_rna"
        )

    # Running PCA for the ADT assay
    # ADT
    markers.oi = intersect(adt.markers, rownames(multimodal[['ADT']]))
    missing = setdiff(adt.markers, markers.oi)
    if (length(missing) > 0) log_message(sprintf('ADT markers not in the data (skipped): %s', paste(missing, collapse = ', ')), 'WARN')
    log_message(sprintf('ADT PCA on %d markers', length(markers.oi)), 'INFO')

    DefaultAssay(multimodal) = 'ADT'
    multimodal = JoinLayers(multimodal, layers = c('counts', 'data', 'scale.data'))

    DefaultAssay(multimodal) = 'ADT'
    multimodal = multimodal %>%
        ScaleData() %>% 
        RunPCA(
            assay = "ADT",
            npcs = adt.npcs,
            reduction.key = "adtPC_", 
            reduction.name = "pca_adt",
            features = markers.oi
        )

    # -----------------------------------------------------------------------------------------------

    return(multimodal)
}

#' Harmony on the RNA PCA (reduction harmony_rna)
#'
#' @param merged.seurat Seurat object with pca_rna.
#' @param group.by.vars Character. Batch variables (orig.ident before).
#'
#' @return The Seurat object with harmony_rna.
harmonize_rna_assay = function(merged.seurat, group.by.vars) {

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
    set.seed(seed)
    harmonized.seurat = RunHarmony(
        merged.seurat,
        group.by.vars = group.by.vars, 
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

#' Harmony on the ADT PCA (reduction harmony_adt)
#'
#' @param merged.seurat Seurat object with pca_adt.
#' @param group.by.vars Character. Batch variables.
#'
#' @return The Seurat object with harmony_adt.
harmonize_adt_assay = function(merged.seurat, group.by.vars) {
  
  # Manually set variable features of merged Seurat object
  DefaultAssay(merged.seurat) = 'ADT'
  
  # Running Harmony
  # -----------------------------------------------------------------------------------------------
  # Run Harmony to integrate samples
  # By default, harmony will look up the pca cell embeddings and use these to run harmony.
  # Therefore, it assumes that the Seurat object has these embeddings already precomputed.
  # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
    set.seed(seed)
    harmonized.seurat = RunHarmony(
        merged.seurat,
        group.by.vars = group.by.vars, 
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

# -- MAIN ----------------------------------------------------------------------
seu <- qs2::qs_read(file = smk@input$merged_seurat)

adt_markers <- if (!is.null(h$adt_markers_file) && nzchar(h$adt_markers_file)) {
  read.csv(h$adt_markers_file, header = FALSE)$V1
} else {
  grep(smk@params$isotype_pattern, rownames(seu[["ADT"]]), value = TRUE, invert = TRUE)
}

seu <- process_multimodal_seurat(seu, adt.markers = adt_markers,
                                 rna.npcs = as.integer(h$rna_npcs), adt.npcs = as.integer(h$adt_npcs))
seu <- harmonize_rna_assay(seu, group.by.vars = unlist(h$group_by_vars))
seu <- harmonize_adt_assay(seu, group.by.vars = unlist(h$group_by_vars))

DefaultAssay(seu) <- "RNA"
seu <- run_wnn(seu, k.nn.val = wnn$k_nn, prune.SNN.val = wnn$prune_snn,
               rna.dims = wnn$rna_dims, adt.dims = wnn$adt_dims)
seu <- run_wnn_umap(seu)

qs2::qs_save(object = seu, file = smk@output$harmonized_seurat)
log_message("Integration done", "DONE")
