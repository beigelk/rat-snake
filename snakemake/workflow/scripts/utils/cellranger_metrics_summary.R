library(tidyverse)

read_metrics <- function(f) {
    sample <- str_remove(basename(f), "__metrics_summary.csv")
    read_csv(f, show_col_types = FALSE) |> mutate(sample = sample)
}

metrics <- list.files("metrics_summaries",
                       pattern = "*.csv",
                       full.names = TRUE) |>
    map_dfr(read_metrics)

summary_tsv <- metrics |>
    mutate(
        col_name = case_when(
            # For rows with a meaningful Group Name, collapse all four fields
            !is.na(`Group Name`) & `Grouped By` != "" ~
                paste(`Category`, `Library Type`, `Grouped By`, `Metric Name`, sep = " - "),
            # For Cells rows with no grouping
            TRUE ~
                paste(`Category`, `Library Type`, `Metric Name`, sep = " - ")
        )
    ) |>
    select(sample, col_name, `Metric Value`) |>
    # If there are duplicate col_name per sample (multiple fastq IDs etc), take first
    group_by(sample, col_name) |>
    slice(1) |>
    ungroup() |>
    pivot_wider(names_from = col_name, values_from = `Metric Value`)

write_tsv(summary_tsv, "cellranger_metrics_summary.tsv")
