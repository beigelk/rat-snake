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
    paste0('snakemake/workflow/scripts/utils/', c('file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)


# -- FUNCTIONS ------------------------------ #
run_wnn = function(
    seurat.obj,
    snn.graph.name = "wsnn_harmony",
    rna.dims = 15,
    adt.dims = 15
) {

    if (snn.graph.name == "wsnn_harmony") {
        knn.graph.name = "wknn_harmony"
        weighted.nn.name = "weighted.nn_harmony"
    } else {
        message(
            paste0("Unknown graph name. Choose from the following options: ",
            paste0(c(knn.graph.name), sep = ", "))
        )
    }

    # Run WNN on Harmony reduction
    seurat.obj = FindMultiModalNeighbors(
        object = seurat.obj,
        reduction.list = list("harmony_rna", "harmony_adt"),
        knn.graph.name = knn.graph.name,
        snn.graph.name = snn.graph.name,
        weighted.nn.name = weighted.nn.name,
        dims.list = list(1:rna.dims, 1:adt.dims),
        modality.weight.name = c("RNA.weight", 'ADT.weight'),
        k.nn = 10,
        prune.SNN = 1/30,
        verbose = TRUE
    )
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
seu = qs2::qs_read(
    file = seurat_obj
)

# set assay
DefaultAssay(seu) = 'RNA'

# -- READ IN CLUSTER RESULTS ---------------- #
clusters = readr::read_tsv(cluster_results)


# -- ADD CLUSTERS TO OBJECT ----------------- #
if (all.equal(colnames(seu), clusters[['cell_barcode']])) {
    seu = AddMetaData(
        seu,
        metadata = clusters
    )
} else {
    message(
        paste0(
            "Number of cells in Seurat object: ", length(colnames(seu)), "\n",
            "Number of cells in cluster annotation dataframe: ", nrow(clusters)
        )
    )
    stop()
}


# -- SUBSET OBJECT -------------------------- #
message(paste0("Number of cells before subclustering: ", ncol(seu)))
cluster_col = paste0('cluster__', base_name, '__', base_graph_name, '_res.', base_resolution)
message(paste0("Subsetting to keep clusters ", paste0(base_clusters, sep = ' '), paste0(' in cluster scheme ',  cluster_col, '\n')))
seu = subset(
    seu,
    subset = !!sym(cluster_col) %in% base_clusters
)
message(paste0("Number of cells after subclustering: ", ncol(seu)))


# -- REANALYZE ----------------------------------
# Find multimodal neighbors and calculate wnn
print('wnn')
seu = run_wnn(
    seurat.obj = seu,
    snn.graph.name = graph_name
)

# Calculate UMAP from WNN
print('umap')
seu = run_umap(seu)


# -- SAVE OUTPUT ---------------------------- #
qs2::qs_save(
    seu,
    file = seurat_sub_obj
)

# -- SESSION INFO ---------------------------- #
sessionInfo()
