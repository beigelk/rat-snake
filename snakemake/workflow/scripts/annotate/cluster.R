# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj      = snakemake@input[['seurat_obj']]

# params
run_config       = snakemake@params[['run_config']]
pop_name         = snakemake@params[['pop_name']]
graph_name       = snakemake@params[['graph_name']]
resolution       = snakemake@params[['resolution']]
r_library        = snakemake@params[['r_library']]

# output
results_tsv = snakemake@output[['results_tsv']]

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

# utility
library(qs2)
library(tidyverse)
library(readr)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 200000 * 1024^2)
# plan("multicore", workers = 16)


# -- SOURCE SCRIPTS ------------------------- #
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)

# -- FUNCTIONS ------------------------------ #
run_wnn = function(seurat.obj) {
    rna_dims = 15
    adt_dims = 15
    # Run WNN on Harmony reduction
    seurat.obj = FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_adt"),
        knn.graph.name = "wknn_harmony",
        snn.graph.name = "wsnn_harmony",
        weighted.nn.name = "weighted.nn_harmony",
        dims.list = list(1:rna_dims, 1:adt_dims),
        modality.weight.name = c("RNA.weight", 'ADT.weight'),
        k.nn = 10,
        prune.SNN = 1/30,
        verbose = TRUE
    )
}


run_umap = function(seurat.obj) {
    # Run UMAP for the weighted nearest neighbors
    seurat.obj = RunUMAP(
        seurat.obj,
        nn.name = "weighted.nn_harmony",
        reduction.name = "wnn.umap_harmony",
        reduction.key = "wnnUMAP_"
    )

    return(seurat.obj)
}


# == MAIN =================================== #
# -- LOAD OBJECT ---------------------------- #
seu = qs2::qs_read(
    file = seurat_obj
)

# set assay
DefaultAssay(seu) = 'RNA'


# Find multimodal neighbors and calculate wnn
print('wnn')
seu = run_wnn(seu)

# Calculate UMAP from WNN
print('umap')
seu = run_umap(seu)

# -- FIND CLUSTERS -------------------------- #
# clear old clusters
for (col in colnames(seu@meta.data)[grepl(paste0(graph_name, '_res*'), colnames(seu@meta.data))]) {
    seu@meta.data[, col] = NULL
}

# seu = FindNeighbors(seu, reduction = 'pca_adt')

# # algorithm = 4 is Leiden clustering
# seu = FindClusters(
#     seu,
#     algorithm = 1,
#     graph.name = 'wsnn_harmony',
#     resolution = 0.2
# )

# Find clusters
seu = FindClusters(
    seu,
    algorithm = 1,
    graph.name = graph_name,
    resolution = resolution
)

cluster_col = paste0('cluster__', pop_name, '__', graph_name, '_res.', resolution)
cluster_results <- tibble(
    cell_barcode     = colnames(seu),
    !!(cluster_col) := as.character(seu@meta.data[[paste0(graph_name, '_res.', resolution)]])
)


# -- SAVE OUTPUT ---------------------------- #
readr::write_tsv(
    cluster_results,
    file = results_tsv
)


# -- SESSION INFO ---------------------------- #
sessionInfo()
