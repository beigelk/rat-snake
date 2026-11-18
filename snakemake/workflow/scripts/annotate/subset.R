# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj      = snakemake@input[['seurat_obj']]
cluster_results = snakemake@input[['cluster_results']]

# params
run_config      = snakemake@params[['run_config']]
name            = snakemake@params[['pop_name']]
graph_name      = snakemake@params[['graph_name']]
base_name       = snakemake@params[['base_name']]
base_graph_name = snakemake@params[['base_graph_name']]
base_resolution = snakemake@params[['base_resolution']]
base_clusters   = snakemake@params[['base_clusters']]
r_library       = snakemake@params[['r_library']]

# output
seurat_sub_obj = snakemake@output[['seurat_sub_obj']]

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
options(future.globals.maxSize = 100000 * 1024^2)
# plan("multicore", workers = 16)


# -- SOURCE SCRIPTS ------------------------- #
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('utils.R', 'file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)


# -- FUNCTIONS ------------------------------ #
run_wnn = function(
    seurat.obj,
    snn.graph.name = "wsnn_harmony",
    rna.dims = 15,
    adt.dims = 15,
    k.nn.val = 10,
    prune.SNN.val = 1/30
) {

    if (snn.graph.name == "wsnn_harmony") {
        knn.graph.name = "wknn_harmony"
        weighted.nn.name = "weighted.nn_harmony"
    } else {
        log_msg(
            paste0("Unknown graph name. Choose from the following options: ",
            paste0(c(snn.graph.name), sep = ", "))
        )
        stop()
    }

    log_msg(sprintf("RNA dims: 1:%s",    rna.dims), level =  'INFO')
    log_msg(sprintf("ADT dims: 1:%s",    adt.dims), level =  'INFO')
    log_msg(sprintf("kNN value: %s",     k.nn.val), level =  'INFO')
    log_msg(sprintf("prune SNN val: %s", prune.SNN.val), level =  'INFO')

    # Run WNN on Harmony reduction
    log_msg('FindMultiModalNeighbors starting.')
    seurat.obj = FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_adt"),
        knn.graph.name = knn.graph.name,
        snn.graph.name = snn.graph.name,
        weighted.nn.name = weighted.nn.name,
        dims.list = list(1:rna.dims, 1:adt.dims),
        modality.weight.name = c("RNA.weight", 'ADT.weight'),
        k.nn = k.nn.val,
        prune.SNN = prune.SNN.val,
        verbose = TRUE
    )
    log_msg('FindMultiModalNeighbors complete.')

    return(seurat.obj)
}

run_umap = function(
    seurat.obj,
    nn.name = "weighted.nn_harmony",
    assay = "RNA"
) {
    if (nn.name == "weighted.nn_harmony") {
        reduction.name = "wnn.umap_harmony"
        reduction.key = "wnnUMAP_"
    } else {
        message(
            paste0("Unknown nearset neighbor name. Choose from the following options: ",
            paste0(c(nn.name), sep = ", "))
        )
    }

    seurat.obj = RunUMAP(
        seurat.obj,
        nn.name = nn.name,
        reduction.name = reduction.name,
        reduction.key = reduction.key 
    )
    # ----------------------------------------------------------------
    return(seurat.obj)
}


# == MAIN =================================== #
# -- LOAD OBJECT ---------------------------- #
log_msg("Reading Seurat.")
seu = qs2::qs_read(
    file = seurat_obj
)
log_msg("Seurat loaded.")

# set assay
log_msg("Setting DefaultAssay to RNA.")
DefaultAssay(seu) = 'RNA'

# -- READ IN CLUSTER RESULTS ---------------- #
log_msg("Reading in cluster results for base object to subset from.")
clusters = readr::read_tsv(cluster_results)

# -- ADD CLUSTERS TO OBJECT ----------------- #
if (all.equal(colnames(seu), clusters[['cell_barcode']])) {
    seu = AddMetaData(
        seu,
        metadata = clusters
    )
} else {
    log_msg(
        paste0(
            "Number of cells in Seurat object: ", length(colnames(seu)), "\n",
            "Number of cells in cluster annotation dataframe: ", nrow(clusters)
        ),
        level = 'INFO'
    )
    stop()
}

# -- SUBSET OBJECT -------------------------- #
log_msg(sprintf("Number of cells before subclustering: %s", ncol(seu)), level = 'INFO')
cluster_col = paste0('cluster__', base_name, '__', base_graph_name, '_res.', base_resolution)
if (identical(base_clusters, "all")) {
    log_msg(
        paste0("base_clusters == 'all' — keeping every cluster from scheme ", cluster_col, "; skipping filtering."),
        level = 'INFO'
    )
    # seu unchanged — no subsetting
} else {
    log_msg(
        paste0(
            "Subsetting to keep clusters ",
            paste0(base_clusters, sep = ' '),
            paste0(' in cluster scheme ',  cluster_col, '\n')
        ),
        level = 'INFO'
    )
    seu = subset(
        seu,
        subset = !!sym(cluster_col) %in% base_clusters
    )
}
log_msg(sprintf("Number of cells after subclustering: ", ncol(seu)), level = 'INFO')


# -- REANALYZE ----------------------------------
# Find multimodal neighbors and calculate wnn
log_msg("Running FindMultiModalNeighbors.")
seu = run_wnn(
    seurat.obj = seu
)

# Calculate UMAP from WNN
log_msg("RunUMAP for wnn graph.")
seu = run_umap(seu)


# -- SAVE OUTPUT ---------------------------- #
log_msg("Saving Seurat object.")
qs2::qs_save(
    seu,
    file = seurat_sub_obj
)
log_msg("Seurat object saved.")

# -- SESSION INFO ---------------------------- #
sessionInfo()
