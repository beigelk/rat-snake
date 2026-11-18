###############################################################################
# File:     marker_features.R
# Purpose:  Read the features of one plot marker set (dotplot, markerplots)
###############################################################################

#' Features of a marker set
#'
#' @param spec List. The marker set: `type: topn` (with `top_n`) reads the top
#'   markers of each cluster from a filtered marker table, in the table's
#'   order; `type: genelist` reads one feature per line.
#' @param features_file Character. The filtered marker table or the gene list.
#'
#' @return Character vector of unique features (may be empty).
resolve_features <- function(spec, features_file) {
  if (spec$type == "topn") {
    markers <- readr::read_tsv(features_file, show_col_types = FALSE)
    if (nrow(markers) == 0) return(character())
    top <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = spec$top_n %||% 5) %>%
      dplyr::ungroup()
    return(unique(top$gene))
  }
  if (spec$type == "genelist") {
    genes <- trimws(readLines(features_file))
    return(unique(genes[nzchar(genes)]))
  }
  stop(sprintf("Unknown marker_set type: %s", spec$type))
}

#' Keep the features found in an assay, logging the ones that are not
#'
#' @param seu Seurat object.
#' @param features Character vector.
#' @param assay Character. Assay to look in.
#'
#' @return The features present in the assay.
present_features <- function(seu, features, assay) {
  found <- intersect(features, rownames(seu[[assay]]))
  missing <- setdiff(features, found)
  if (length(missing) > 0) {
    log_message(sprintf("%d feature(s) not in %s: %s", length(missing), assay, paste(missing, collapse = ", ")), "WARN")
  }
  found
}

#' A placeholder figure, written when a marker set has no features to plot
#'
#' Keeps the rule's declared output so the run (and the report) go on.
#'
#' @param filename Character. Output path without extension.
#' @param title,reason Character. Shown on the figure.
placeholder_plot <- function(filename, title, reason) {
  p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = title, subtitle = reason)
  plot_pdf_and_png(p, filename, width = 6, height = 2, png_only = TRUE)
}
