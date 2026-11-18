###############################################################################
# Rule:     cellranger_metrics
# Purpose:  One table of Cell Ranger multi metrics, one row per sample
#
# Description:
#   From utils/cellranger_metrics_summary.R, which read every CSV in a
#   hard-coded "metrics_summaries" folder. It now reads each sample's
#   metrics_summary.csv from its Cell Ranger run (made here or finished
#   earlier) and labels rows with the sample name instead of the file name.
#   Each metric becomes a column named "Category - Library Type - [Grouped By -]
#   Metric Name".
#
# Inputs:   snakemake@params$metrics  named list: sample -> metrics_summary.csv
# Outputs:  snakemake@output$tsv      wide table (sample x metric)
#           snakemake@output$long     long table (sample, metric, value)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages(library(tidyverse))


#' Read one sample's metrics_summary.csv
#'
#' @param path Character. Path to metrics_summary.csv.
#' @param sample Character. Sample name for the rows.
#'
#' @return Tibble of the CSV with a `sample` column.
read_metrics <- function(path, sample) {
  read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
    mutate(sample = sample)
}


#' Name each metric and keep sample, metric and value
#'
#' @param metrics Tibble from read_metrics() for one or more samples.
#'
#' @return Tibble with `sample`, `metric`, `value` (first value when a metric repeats).
metrics_long <- function(metrics) {
  metrics |>
    mutate(
      metric = case_when(
        !is.na(`Group Name`) & !is.na(`Grouped By`) & `Grouped By` != "" ~
          paste(Category, `Library Type`, `Grouped By`, `Metric Name`, sep = " - "),
        TRUE ~ paste(Category, `Library Type`, `Metric Name`, sep = " - ")
      )
    ) |>
    select(sample, metric, value = `Metric Value`) |>
    group_by(sample, metric) |>
    slice(1) |>
    ungroup()
}


files <- unlist(smk@params$metrics)
log_message(sprintf("Reading %d metrics files", length(files)), "INFO")
long <- map2_dfr(files, names(files), read_metrics) |> metrics_long()
wide <- long |> pivot_wider(names_from = metric, values_from = value)

write_tsv(long, smk@output$long)
write_tsv(wide, smk@output$tsv)
log_message(sprintf("Wrote %d samples x %d metrics to %s", nrow(wide), ncol(wide) - 1, smk@output$tsv), "DONE")
