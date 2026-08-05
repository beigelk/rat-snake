# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
raw_tsv       = snakemake@input[['raw_tsv']]

# params
lfc_threshold  = snakemake@params[['lfc_threshold']]
padj_threshold = snakemake@params[['padj_threshold']]
pct1_threshold = snakemake@params[['pct1_threshold']]
prefix         = snakemake@params[['prefix']]
r_library      = snakemake@params[['r_library']]

# output
filtered_tsv  = snakemake@output[['filtered_tsv']]

# log
logfile = snakemake@log[1][[1]]


# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)


# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")


# -- LOAD PACKAGES ------------------------- #
# utility
library(tidyverse)
library(readr)
library(openxlsx)

# -- MAIN ----------------------------------- #
# =========================================== #
markers <- readr::read_tsv(raw_tsv, show_col_types = FALSE)

if (is.null(markers) || nrow(markers) == 0) {
  message("No markers found — writing empty table.")
  filtered <- markers[0, , drop = FALSE] %||% data.frame()
} else {
  filtered <- markers %>%
  dplyr::filter(
    avg_log2FC > lfc_threshold,
    p_val_adj  < padj_threshold,
    pct.1      > pct1_threshold
  )
}


dir.create(dirname(filtered_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(filtered, filtered_tsv, row.names = FALSE, quote = FALSE, sep = '\t')
openxlsx::write.xlsx(filtered, sub("\\.tsv$", ".xlsx", filtered_tsv), rowNames = FALSE)

message(sprintf("Filtered: %d markers retained (from %d raw)", nrow(filtered), nrow(markers)))