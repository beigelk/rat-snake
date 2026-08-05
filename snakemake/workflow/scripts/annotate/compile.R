# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
target_obj       = snakemake@input[['target_obj']]
metadata_obj     = snakemake@input[['metadata_obj']]
cluster_files    = snakemake@input[['cluster_results']]

# params
compile_target    = snakemake@params[['compile_target']]
cluster_meta_json = snakemake@params[['cluster_meta_json']]
r_library         = snakemake@params[['r_library']]

# output
compiled_obj  = snakemake@output[['compiled_obj']]

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
library(jsonlite)


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




# == MAIN =================================== #
# -- PARSE JSON ------------------------------ #
cluster_meta = fromJSON(cluster_meta_json)
drop_cols <- paste0(
    'cluster__', cluster_meta$run_name,
    '__', cluster_meta$graph_name,
    '_res.', cluster_meta$resolution
)

# -- LOAD OBJECT ---------------------------- #
seu = qs2::qs_read(
    file = target_obj
)

# clear old clusters
for (col in drop_cols) {
    seu@meta.data[, col] = NULL
}

# -- GET FILES TO ADD CLUSTERS TO SEURAT ----- #
cluster_files    <- snakemake@input$cluster_results

all_clusters <- NULL
if (!is.null(cluster_files)) {
    cluster_tables    <- lapply(cluster_files, read.delim, sep = "\t")
    names(cluster_tables)    <- sapply(cluster_tables, colnames)[2,]

    all_clusters <- Reduce(
        function(x, y) merge(x, y, by = "cell_barcode", all = TRUE),
        cluster_tables
    )
}

# -- ADD CLUSTERS TO SEURAT ----------------- #
all_clusters = all_clusters %>%
    column_to_rownames('cell_barcode')

all_clusters = all_clusters[colnames(seu), ]

# -- ADD METADATA TO SEURAT ----------------- #
seu = AddMetaData(
    seu,
    metadata = all_clusters
)

# -- FACTOR CLUSTER COLUMNS IN NUMERICAL ORDER ----------------- #
cluster_cols <- colnames(all_clusters)

for (col in cluster_cols) {
    vals <- seu@meta.data[[col]]
    ordered_levels <- as.character(sort(as.numeric(unique(as.character(vals)))))
    seu@meta.data[[col]] <- factor(as.character(vals), levels = ordered_levels)
}

# -- SAVE OUTPUT ---------------------------- #
qs2::qs_save(
    seu,
    file = compiled_obj
)
