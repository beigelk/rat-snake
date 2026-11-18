###############################################################################
# Rule:     filtermarkers
# Purpose:  Keep the markers that pass one threshold set
#
# Description:
#   Unchanged: avg_log2FC > lfc_threshold, p_val_adj < padj_threshold and
#   pct.1 > pct1_threshold. Row order is FindAllMarkers' (by cluster, then
#   p value), which is the order the "top N" marker sets take markers in.
#
# Inputs:   raw_tsv
# Outputs:  filtered_tsv (and an .xlsx copy)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

t <- smk@params$thresholds
markers <- readr::read_tsv(smk@input$raw_tsv, show_col_types = FALSE,
                           col_types = cols(cluster = col_character(), gene = col_character()))

filtered <- markers %>%
  dplyr::filter(
    avg_log2FC > t$lfc_threshold,
    p_val_adj  < t$padj_threshold,
    pct.1      > t$pct1_threshold
  )

write.table(filtered, smk@output$filtered_tsv, row.names = FALSE, quote = FALSE, sep = "\t")
openxlsx::write.xlsx(filtered, sub("\\.tsv$", ".xlsx", smk@output$filtered_tsv), rowNames = FALSE)
log_message(sprintf("%d of %d markers pass", nrow(filtered), nrow(markers)), "DONE")
