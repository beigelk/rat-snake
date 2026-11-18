###############################################################################
# Rule:     metadata_plots
# Purpose:  UMAPs and QC feature plots of the integrated object, coloured by
#           sample metadata (annotate.metadata.plots)
#
# Description:
#   The plotting half of visualization/metadata.R (its metadata half is the
#   annotate module's metadata rule). Plot functions and the config format
#   (umap / feature_stat, each with defaults and a list of plots) are
#   unchanged; the cell metadata comes from the metadata rule's table.
#   Fix: the script saved its metadata table as qs2 and then overwrote the
#   same file with saveRDS.
#
# Inputs:   seurat_obj (integrated), metadata_obj (cell metadata table)
# Outputs:  figures (folder), index (TSV: plot, file) for the report
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(qs2)
  library(tidyverse)
  library(ggplot2)
})

plot_config <- smk@params$plot_config
figdir <- smk@output$figures
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
figbase <- file.path(figdir, paste0(smk@params$project, "__"))

# -- FUNCTIONS -----------------------------------------------------------------
#' Call a function with only the arguments it accepts
#'
#' @param fn Function.
#' @param args Named list of arguments (extra ones are dropped).
#'
#' @return Whatever `fn` returns.
safe_call <- function(fn, args) {
    valid_args <- names(formals(fn))
    do.call(fn, args[names(args) %in% valid_args])
}

#' UMAP of the cells coloured by one metadata column (scplotter::CellDimPlot)
#'
#' @param seu Seurat object.
#' @param group_var Character. Column to colour by.
#' @param split_var Character or NULL. Column to split panels by.
#' @param subset_var,subset_val Character or NULL. Plot only cells where subset_var == subset_val.
#' @param reduction Character. Reduction to plot.
#' @param width,height Numeric. Size in inches.
#' @param filename_base Character. Start of the output path.
#'
#' @return The plot, invisibly. Writes UMAP__{suffix}.pdf/.png.
plot_umap <- function(
    seu,
    group_var,
    split_var = NULL,
    subset_var = NULL,
    subset_val = NULL,
    reduction = "wnn.umap_harmony",
    width = 12,
    height = 8,
    filename_base
) {
    # safely coerce to NULL
    null_if_missing <- function(x) if (is.null(x) || length(x) == 0 || identical(x, "null") || is.na(x)) NULL else x

    split_var  <- null_if_missing(split_var)
    subset_var <- null_if_missing(subset_var)
    subset_val <- null_if_missing(subset_val)

    obj <- if (!is.null(subset_var)) {
        subset(seu, subset = !!sym(subset_var) == subset_val)
    } else {
        seu
    }

    p <- if (is.null(split_var) || group_var == split_var) {
        CellDimPlot(obj, group_by = group_var, reduction = reduction)
    } else {
        CellDimPlot(obj, group_by = group_var, split_by = split_var, reduction = reduction)
    }

    suffix <- paste(
        c(
            subset_val,
            if (is.null(split_var) || group_var == split_var) group_var
            else paste0(group_var, "__", split_var)
        ),
        collapse = "_"
    )

    plot_pdf_and_png(
        p,
        paste0(filename_base, "UMAP__", suffix),
        width  = width,
        height = height
    )
    invisible(p)
}

#' Violin or ridge plot of QC features by one metadata column (scplotter::FeatureStatPlot)
#'
#' @param seu Seurat object.
#' @param features Character. Features (metadata columns or genes) to plot, one plot each.
#' @param group_var Character. Column on the x axis.
#' @param split_var Character or NULL. Unused by the current plots; kept for the config.
#' @param plot_type Character. "violin" or "ridge".
#' @param subset_var,subset_val Character or NULL. Plot only cells where subset_var == subset_val.
#' @param width,height Numeric. Size in inches.
#' @param filename_base Character. Start of the output path.
#'
#' @return NULL. Writes {VlnPlot|RidgePlot}__{suffix}.pdf/.png per feature.
plot_feature_stat <- function(
    seu,
    features,
    group_var,
    split_var = NULL,
    plot_type = "violin",
    subset_var = NULL,
    subset_val = NULL,
    width = 20,
    height = 8,
    filename_base
) {
    # safely coerce to NULL
    null_if_missing <- function(x) if (is.null(x) || length(x) == 0 || identical(x, "null") || is.na(x)) NULL else x

    split_var  <- null_if_missing(split_var)
    subset_var <- null_if_missing(subset_var)
    subset_val <- null_if_missing(subset_val)

    obj <- if (!is.null(subset_var)) {
        subset(seu, subset = !!sym(subset_var) == subset_val)
    } else {
        seu
    }

    prefix <- switch(plot_type, violin = "VlnPlot", ridge = "RidgePlot", plot_type)

    for (feature in features) {
        p <- FeatureStatPlot(
            obj,
            plot_type = plot_type,
            group_by  = group_var,
            features  = feature,
            add_point = TRUE
        ) + NoLegend()

        suffix <- paste(c(subset_val, feature, group_var), collapse = "_")
        plot_pdf_and_png(
            p,
            paste0(filename_base, prefix, "__", suffix),
            width  = width,
            height = height
        )
    }
}


# -- MAIN ----------------------------------------------------------------------
seu <- qs2::qs_read(file = smk@input$seurat_obj)
meta <- qs2::qs_read(smk@input$metadata_obj)
new_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
seu <- AddMetaData(seu, metadata = meta[colnames(seu), new_cols, drop = FALSE])

umap_cfg <- plot_config$umap
if (!is.null(umap_cfg)) {
  defaults <- modifyList(list(reduction = smk@params$reduction), umap_cfg$defaults %||% list())
  for (p in umap_cfg$plots) {
    args <- modifyList(defaults, p)
    safe_call(plot_umap, c(list(seu = seu, filename_base = figbase), args))
  }
}

fs_cfg <- plot_config$feature_stat
if (!is.null(fs_cfg)) {
  features <- intersect(unlist(fs_cfg$features), colnames(seu@meta.data))
  for (p in fs_cfg$plots) {
    args <- modifyList(fs_cfg$defaults %||% list(), p)
    safe_call(plot_feature_stat, c(list(seu = seu, features = features, filename_base = figbase), args))
  }
}

pngs <- list.files(figdir, pattern = "\\.png$", full.names = TRUE)
write_tsv(tibble(plot = tools::file_path_sans_ext(basename(pngs)), file = pngs), smk@output$index)
log_message(sprintf("%d metadata plots", length(pngs)), "DONE")
