
# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj      = snakemake@input[['seurat_obj']]

# params
metadata_json    = snakemake@params[['metadata_json']]
plot_config      = snakemake@params[['plot_config']]
figbase          = snakemake@params[['figbase']]
tblbase          = snakemake@params[['tblbase']]
r_library        = snakemake@params[['r_library']]

# output
metadata_obj = snakemake@output[['metadata_obj']]

# log
logfile = snakemake@log[1][[1]]


# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)


# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")


# -- VIEW INPUTS --------------------------- #
dput(snakemake)


# -- LOAD PACKAGES ------------------------- #
# analysis
library(Seurat)
library(scplotter)
library(scCustomize)

# utility
library(qs2)
library(tidyverse)
library(readxl)
library(stringdist)
library(jsonlite)

# plot
library(ggplot2)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
# options(future.globals.maxSize = 200000 * 1024^2)


# -- SOURCE SCRIPTS ------------------------- #
sapply(c('snakemake/workflow/scripts/utils/plot.R'), source, print.eval = FALSE)


# -- FUNCTIONS ------------------------------ #
# -- Safe caller: only passes args the function accepts
safe_call <- function(fn, args) {
    valid_args <- names(formals(fn))
    do.call(fn, args[names(args) %in% valid_args])
}

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


# -- MAKE OUPUT DIRS ------------------------ #
figdir = dirname(figbase)
tbldir = dirname(tblbase)
objdir = dirname(metadata_obj)
sapply(
    c(figdir, tbldir, objdir),
    dir.create, showWarnings = FALSE, recursive = TRUE
)


# == MAIN =================================== #
# -- LOAD OBJECT ---------------------------- #
seu = qs2::qs_read(
    file = seurat_obj
)

# Clean names and add as new metadata in Seurat object
seu@meta.data$Subject_ID = gsub("[ -]+", "_", seu@meta.data$sample.nameID)

# -- LOAD METADATA -------------------------- #
# metadata is now a named list, one entry per dataset
# metadata$dataset1$file, metadata$dataset1$cols, etc.
metadata_config <- fromJSON(metadata_json)

# -- ADD METADATA --------------------------- #
metadata_list = vector(mode = 'list')
for (metadata_input in metadata_config) {

    metadata = readxl::read_excel(
        path = metadata_input$file,
        sheet = metadata_input$sheet,
        n_max = metadata_input$nrows
    )

    metadata = metadata %>%
        dplyr::select(all_of(c(metadata_input$id_col, metadata_input$cols))) %>%
        mutate(Subject_ID = gsub("[ -]+", "_", .data[[metadata_input$id_col]])) %>% 
        mutate(across(everything(), ~ {
            converted = suppressWarnings(as.numeric(.x))
            if (all(is.na(converted) == is.na(.x))) converted else .x
        })) %>%
        mutate(
            Dataset = metadata_input$name
        )

    print(metadata)

    levels = metadata_input$levels
    for (col in names(levels)) {
        if (col %in% colnames(metadata)) {
            lvls = metadata_input$levels[[col]]
            metadata[[col]] = factor(
                metadata[[col]],
                levels = lvls
            )
        } else {
            message(paste0(col, ' in level config not present in metadata table for dataset ', metadata_input$name))
        }
    }

    metadata_list[[basename(metadata_input$name)]] = list(df = metadata, id_col = 'Subject_ID')

}

metadata_joined = bind_rows(lapply(metadata_list, `[[`, "df"))
colnames(metadata_joined)

metadata_joined = metadata_joined %>%
    mutate(Subject_ID = gsub("[ -]+", "_", Subject_ID))

tbl_meta_ids = metadata_joined$Subject_ID
seu_meta_ids = unique(seu@meta.data$Subject_ID)

matches = stringdist::amatch(seu_meta_ids, tbl_meta_ids, method = "jw", maxDist = 0.05)
matches_df = data.frame(
    metadata_table = tbl_meta_ids[matches],
    seurat_metadata = seu_meta_ids,
    distance = stringdist::stringdist(seu_meta_ids, tbl_meta_ids[matches], method = "jw")
)

new_metadata = metadata_joined %>%
    mutate(Original_Subject_ID = Subject_ID) %>%
    mutate(Subject_ID = ifelse(
        Subject_ID %in% matches_df$metadata_table,
        matches_df$seurat_metadata[match(Subject_ID, matches_df$metadata_table)],
        Subject_ID
    )) %>%
    relocate(Subject_ID) %>%
    relocate(Original_Subject_ID, .after = Subject_ID)

seu_meta = seu@meta.data %>%
    left_join(
        new_metadata,
        by = join_by(Subject_ID == Subject_ID)
    ) %>%
    dplyr::select(all_of(c(colnames(new_metadata))))
head(seu_meta)

seu = AddMetaData(
    seu,
    metadata = seu_meta
)


# -- METADATA PLOTS ------------------------- #
features <- c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent_mt", "IsoCount")

# -- UMAPs ---------------------------------------------------------- ##
umap_cfg  <- plot_config$umap
defaults  <- umap_cfg$defaults

for (p in umap_cfg$plots) {
    args <- modifyList(defaults, p)
    safe_call(plot_umap, c(list(seu = seu, filename_base = figbase), args))
}

# -- Feature stat plots --------------------------------------------- ##
fs_cfg   <- plot_config$feature_stat
defaults <- fs_cfg$defaults
features <- fs_cfg$features

for (p in fs_cfg$plots) {
    args <- modifyList(defaults, p)
    
    safe_call(
        plot_feature_stat,
        c(list(seu = seu, features = features, filename_base = figbase), args))
}


# -- SAVE OUTPUT ---------------------------- #
qs2::qs_save(
    seu@meta.data,
    file = metadata_obj
)

saveRDS(
    seu@meta.data,
    file = metadata_obj
)
