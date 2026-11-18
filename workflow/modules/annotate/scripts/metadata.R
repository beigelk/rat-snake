###############################################################################
# Rule:     metadata
# Purpose:  Build the cell metadata table: sample metadata from each dataset's
#           spreadsheet and, with TCR data, scRepertoire clonotypes
#
# Description:
#   Sample metadata (unchanged logic): each dataset's spreadsheet is read
#   (id_col + cols, optional factor levels), sample names are cleaned
#   ("[ -]+" -> "_") into Subject_ID and matched to the samples by
#   Jaro-Winkler distance (annotate.metadata.match_max_dist, 0.05 before),
#   and the columns are joined to every cell.
#   Fixes: the joined table kept its cell names (left_join dropped the row
#   names, so AddMetaData filled every column with NA); the plots moved to
#   the visualization module's metadata_plots rule; CSV/TSV sheets work too.
#   Clonotypes: scRepertoire::combineExpression() adds CTgene, CTnt, CTaa,
#   CTstrict, clonalProportion/clonalFrequency and cloneSize to each cell.
#
# Inputs:   seurat_obj (integrated), clones (combineTCR output per sample)
# Outputs:  metadata_obj (qs2 data frame, rows = cells), tsv (same table),
#           matches (TSV: how sample names were matched to the spreadsheets)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(readxl)
  library(stringdist)
})

p <- smk@params


#' Clean a sample name for matching ("[ -]+" -> "_")
#'
#' @param x Character vector.
#'
#' @return Character vector.
clean_id <- function(x) gsub("[ -]+", "_", x)


#' Read one dataset's sample metadata sheet
#'
#' @param m List. The dataset's metadata block (file, sheet, id_col, cols,
#'   nrows, name, levels).
#'
#' @return Tibble with Subject_ID, the chosen columns (numbers converted where
#'   every value is a number), Dataset, and factors for `levels`.
read_sample_metadata <- function(m) {
  ext <- tolower(tools::file_ext(m$file))
  nmax <- if (is.null(m$nrows)) Inf else m$nrows
  tbl <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path = m$file, sheet = m$sheet, n_max = nmax)
  } else {
    readr::read_delim(m$file, delim = if (ext == "csv") "," else "\t", n_max = nmax, show_col_types = FALSE)
  }
  tbl <- tbl %>%
    dplyr::select(all_of(unique(c(m$id_col, unlist(m$cols))))) %>%
    mutate(Subject_ID = clean_id(.data[[m$id_col]])) %>%
    mutate(across(everything(), ~ {
      converted <- suppressWarnings(as.numeric(.x))
      if (all(is.na(converted) == is.na(.x))) converted else .x
    })) %>%
    mutate(Dataset = m$name)
  for (col in names(m$levels)) {
    if (col %in% colnames(tbl)) {
      tbl[[col]] <- factor(tbl[[col]], levels = unlist(m$levels[[col]]))
    } else {
      log_message(sprintf("%s in levels is not a column of the %s metadata", col, m$name), "WARN")
    }
  }
  tbl
}


#' Match the samples' Subject_IDs to the metadata table's IDs
#'
#' @param seu_ids Character. Subject_IDs in the object.
#' @param tbl_ids Character. Subject_IDs in the metadata.
#' @param max_dist Numeric. Largest Jaro-Winkler distance accepted.
#'
#' @return Data frame: seurat_metadata, metadata_table (NA when unmatched), distance.
match_ids <- function(seu_ids, tbl_ids, max_dist) {
  idx <- stringdist::amatch(seu_ids, tbl_ids, method = "jw", maxDist = max_dist)
  data.frame(
    seurat_metadata = seu_ids,
    metadata_table  = tbl_ids[idx],
    distance        = stringdist::stringdist(seu_ids, tbl_ids[idx], method = "jw")
  )
}


# -- MAIN ----------------------------------------------------------------------
seu <- qs2::qs_read(file = smk@input$seurat_obj)
seu$Subject_ID <- clean_id(seu$Sample_ID)
meta <- seu@meta.data
matches_df <- data.frame(seurat_metadata = unique(meta$Subject_ID), metadata_table = NA, distance = NA)

sheets <- Filter(Negate(is.null), p$metadata)
if (length(sheets) > 0) {
  metadata_joined <- bind_rows(lapply(sheets, read_sample_metadata)) %>%
    mutate(Subject_ID = clean_id(Subject_ID))
  matches_df <- match_ids(unique(meta$Subject_ID), metadata_joined$Subject_ID, p$match_max_dist)
  unmatched <- matches_df$seurat_metadata[is.na(matches_df$metadata_table)]
  if (length(unmatched) > 0) {
    log_message(sprintf("No metadata row for: %s", paste(unmatched, collapse = ", ")), "WARN")
  }

  # Rename the metadata's IDs to the object's, then join to every cell
  new_metadata <- metadata_joined %>%
    mutate(Original_Subject_ID = Subject_ID) %>%
    mutate(Subject_ID = ifelse(
      Subject_ID %in% matches_df$metadata_table,
      matches_df$seurat_metadata[match(Subject_ID, matches_df$metadata_table)],
      Subject_ID
    )) %>%
    relocate(Subject_ID) %>%
    relocate(Original_Subject_ID, .after = Subject_ID)

  added <- setdiff(colnames(new_metadata), "Subject_ID")
  joined <- meta %>%
    select(-any_of(setdiff(added, "Dataset"))) %>%
    rownames_to_column("cell") %>%
    left_join(new_metadata, by = "Subject_ID", suffix = c("", ".metadata")) %>%
    column_to_rownames("cell")
  if ("Dataset.metadata" %in% colnames(joined)) {
    joined$Dataset <- dplyr::coalesce(as.character(joined$Dataset.metadata), as.character(joined$Dataset))
    joined$Dataset.metadata <- NULL
  }
  meta <- joined[rownames(meta), , drop = FALSE]
}

if (length(smk@input$clones) > 0) {
  suppressPackageStartupMessages(library(scRepertoire))
  tcr <- p$tcr
  clones <- unlist(lapply(unlist(smk@input$clones), qs2::qs_read), recursive = FALSE)
  log_message(sprintf("Adding clonotypes from %d samples", length(clones)), "INFO")
  seu@meta.data <- meta
  args <- list(input.data = clones, sc.data = seu, cloneCall = tcr$clone_call, proportion = isTRUE(tcr$proportion))
  if (!is.null(tcr$clone_size_bins)) args$cloneSize <- unlist(tcr$clone_size_bins)
  seu <- do.call(combineExpression, args)
  meta <- seu@meta.data
  log_message(sprintf("%d of %d cells have a clonotype", sum(!is.na(meta$CTstrict)), nrow(meta)), "INFO")
}

qs2::qs_save(meta, file = smk@output$metadata_obj)
write_tsv(rownames_to_column(meta, "cell"), smk@output$tsv)
write_tsv(matches_df, smk@output$matches)
log_message(sprintf("Metadata for %d cells, %d columns", nrow(meta), ncol(meta)), "DONE")
