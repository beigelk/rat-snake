# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj = snakemake@input[['seurat_obj']]

# params
pop_name       = snakemake@params[['pop_name']]
graph_name     = snakemake@params[['graph_name']]
assay_name     = snakemake@params[['assay']]
resolution     = snakemake@params[['resolution']]
resolution_str = snakemake@params[['resolution_str']]
outfile_base   = snakemake@params[['outfile_base']]
lfc_threshold  = snakemake@params[['lfc_threshold']]
padj_threshold = snakemake@params[['padj_threshold']]
pct1_threshold = snakemake@params[['pct1_threshold']]
prefix         = snakemake@params[['prefix']]
r_library      = snakemake@params[['r_library']]

# output
donefile = snakemake@output[['done']]

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
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# Error: package or namespace load failed for ‘AnnotationDbi’ in dyn.load(file, DLLpath = DLLpath, ...):
#  unable to load shared object 'renv/library/linux-rhel-9.8/R-4.5/x86_64-pc-linux-gnu/RSQLite/libs/RSQLite.so':
#   openldap/2.6.4/lib/libldap.so.2: undefined symbol: EVP_md2, version OPENSSL_3.0.0
# Error: package ‘AnnotationDbi’ could not be loaded
# Execution halted

# utility
library(qs2)
library(tidyverse)
library(readr)
library(openxlsx)

# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 200000 * 1024^2)
# plan("multicore", workers = 16)


# -- SOURCE SCRIPTS ------------------------- #
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('file_io.R', 'utils.R')),
    source,
    print.eval = FALSE
)

# -- FUNCTIONS ------------------------------ #
find_all_markers = function(seu_obj, meta_col, assay_name) {
  markers = tryCatch(
    {
      FindAllMarkers(
        seu_obj,
        assay = assay_name,
        group.by = meta_col,
        only.pos = TRUE,
        verbose = TRUE
      )
    },
    error = function(e) {
      message('FindMarkers failed: ', e$message)
      return(NULL)
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    return(NULL)
  }

  markers = markers %>%
    remove_rownames() %>%
    dplyr::select(c(cluster, gene, pct.1, pct.2, avg_log2FC, p_val, p_val_adj))

  # markers$entrez <- AnnotationDbi::mapIds(
  #   org.Hs.eg.db, keys = markers$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  # )

  # entrez_keys <- as.character(markers$entrez)

  # markers$gene_name <- AnnotationDbi::mapIds(
  #   org.Hs.eg.db,
  #   keys = ifelse(is.na(entrez_keys), "1", entrez_keys),
  #   column = "GENENAME", keytype = "ENTREZID", multiVals = "first"
  # ) %>% unname()

  # markers$gene_name[is.na(markers$entrez)] <- NA

  return(markers)
}

# # Alternate method of writing out markes, original for running FindAllMarkers in a loop per res
# write_markers = function(markers, output_dir, prefix,
#                           lfc_threshold = 1, padj_threshold = 0.05, pct1_threshold = 0.1) {
#   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#   if (is.null(markers) || nrow(markers) == 0) {
#     message("No markers to write for: ", prefix)
#     return(invisible(NULL))
#   }

#   raw_path <- file.path(output_dir, paste0(prefix, '__markers_raw'))
#   write.table(markers, paste0(raw_path, '.tsv'), row.names = FALSE, quote = FALSE, sep = '\t')
#   openxlsx::write.xlsx(markers, paste0(raw_path, '.xlsx'), rowNames = FALSE)
#   message("Written: ", raw_path)

#   filtered <- markers %>%
#     dplyr::filter(avg_log2FC > lfc_threshold, p_val_adj < padj_threshold, pct.1 > pct1_threshold)

#   filtered_path <- file.path(output_dir, paste0(prefix, '__markers_filtered'))
#   write.table(filtered, paste0(filtered_path, '.tsv'), row.names = FALSE, quote = FALSE, sep = '\t')
#   openxlsx::write.xlsx(filtered, paste0(filtered_path, '.xlsx'), rowNames = FALSE)
#   message("Written: ", filtered_path, " (", nrow(filtered), " markers retained)")

#   invisible(list(raw = markers, filtered = filtered))
# }

# -- MAIN ------------------------------------------ #
# ================================================== #
seu <- qs2::qs_read(seurat_obj)

cluster_col <- sprintf("cluster__%s__%s_res.%s", pop_name, graph_name, resolution_str)
if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(sprintf("Column '%s' not found in metadata.", cluster_col))
}

Idents(seu) <- cluster_col
DefaultAssay(seu) <- assay_name

markers <- find_all_markers(seu_obj = seu, meta_col = cluster_col, assay_name = assay_name)

raw_tsv <- snakemake@output[['raw_tsv']]
dir.create(dirname(raw_tsv), recursive = TRUE, showWarnings = FALSE)

if (is.null(markers) || nrow(markers) == 0) {
  message("No markers found — writing empty table.")
  markers <- markers[0, , drop = FALSE] %||% data.frame()
}

write.table(markers, raw_tsv, row.names = FALSE, quote = FALSE, sep = '\t')
openxlsx::write.xlsx(markers, sub("\\.tsv$", ".xlsx", raw_tsv), rowNames = FALSE)

# sink(type = "message")
# sink(type = "output")
# ================================================== #

