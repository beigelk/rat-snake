###############################################################################
# Rule:     tcr_clones
# Purpose:  Build one sample's TCR clonotypes with scRepertoire
#
# Description:
#   Wired from the unfinished clones.R. Reads the sample's
#   filtered_contig_annotations.csv from Cell Ranger multi and runs
#   scRepertoire::combineTCR() with the original settings (removeNA,
#   removeMulti and filterMulti all FALSE by default). Barcodes become
#   "{sample}_{barcode}", which matches the cell names after the merge step,
#   so the annotate module's metadata rule can add the clonotypes to every
#   cell with combineExpression().
#
# Inputs:   snakemake@params$contigs   filtered_contig_annotations.csv
# Outputs:  clones      qs2: combineTCR() output (a list with one data frame)
#           table       TSV of the same data frame
#           figure      PNG: cells per clonotype
#           summary     TSV: sample, tcr_cells, clonotypes, largest_clone
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(scRepertoire)
  library(qs2)
  library(tidyverse)
  library(ggplot2)
})

sample_name <- smk@wildcards$sample
opts <- smk@params$opts
clone_col <- paste0("CT", opts$clone_call)


contigs <- read.csv(smk@params$contigs)
log_message(sprintf("%d contigs for %d barcodes", nrow(contigs), length(unique(contigs$barcode))), "INFO")

combined <- combineTCR(
  setNames(list(contigs), sample_name),
  samples     = sample_name,
  removeNA    = isTRUE(opts$remove_na),
  removeMulti = isTRUE(opts$remove_multi),
  filterMulti = isTRUE(opts$filter_multi)
)

clones <- combined[[1]]
qs2::qs_save(combined, file = smk@output$clones)
write_tsv(clones, smk@output$table)

sizes <- clones %>%
  filter(!is.na(.data[[clone_col]])) %>%
  count(.data[[clone_col]], name = "cells")

p <- ggplot(sizes, aes(x = cells)) +
  geom_histogram(binwidth = 1, boundary = 0.5) +
  scale_y_log10() +
  labs(title = sprintf("%s: cells per clonotype (%s)", sample_name, clone_col),
       x = "Cells in clonotype", y = "Clonotypes (log10)") +
  theme_bw()
plot_pdf_and_png(p, tools::file_path_sans_ext(smk@output$figure), width = 7, height = 5)

write_tsv(
  tibble(sample = sample_name, tcr_cells = nrow(clones), clonotypes = nrow(sizes),
         largest_clone = if (nrow(sizes)) max(sizes$cells) else 0L),
  smk@output$summary
)
log_message(sprintf("%d cells with a TCR in %d clonotypes", nrow(clones), nrow(sizes)), "DONE")
