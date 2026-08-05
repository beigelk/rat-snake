# -- SNAKEMAKE OBJ FOR DEBUGGING ----------- #
saveRDS(snakemake, file = snakemake@params[['smk_obj']])

# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj      = snakemake@input[['seurat_obj']]

# params
metadata_json    = snakemake@params[['metadata_json']]
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
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)


# -- FUNCTIONS ------------------------------ #


# -- MAKE OUPUT DIRS ------------------------ #
objdir = dirname(metadata_obj)
sapply(
    c(objdir),
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

# -- SAVE OUTPUT ---------------------------- #
qs2::qs_save(
    seu@meta.data,
    file = metadata_obj
)

saveRDS(
    seu@meta.data,
    file = replace_ext(path = metadata_obj, new_ext = "rds")
)




# p1 = DimPlot(
#         object = seu,
#         group.by = 'Sample_ID',
#         reduction = 'wnn.umap_harmony',
#         shuffle = TRUE) +
#     guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

# plot_pdf_and_png(
#     plot = p1,
#     filename = paste0(filename_base, 'UMAP'),
#     height = 10,
#     width = 12
# )



# # -- UMAP: Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Dataset"
# plot_pdf_and_png(
#     plot = DimPlot(seu, split.by = split_var, group.by = group_var), 
#     filename = paste0(
#         figdir, '/', filename_base, 'UMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 12,
#     height = 8
# )

# # -- UMAP: Sample_ID ---------------------------- ##
# split_var = "Sample_ID"
# group_var = "Sample_ID"
# p = if (group_var == split_var) {
#     DimPlot(seu, group.by = group_var) 
# } else {
#     DimPlot(seu, split.by = split_var, group.by = group_var)
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'UMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 12,
#     height = 8
# )

# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sample_ID"

# p = if (group_var == split_var) {
#     DimPlot(seu, group.by = group_var) 
# } else {
#     DimPlot(seu, split.by = split_var, group.by = group_var)
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'UMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 16,
#     height = 8
# )



# for (sample in unique(seu@meta.data$Sample_ID)) {
#     Idents(seu) = 'Sample_ID'
#     p1 = Cluster_Highlight_Plot(
#         seurat_object = seu,
#         cluster_name = sample,
#         highlight_color = "navy",
#         background_color = "lightgray"
#     )

#     plot_pdf_and_png(
#         plot = p1, 
#         filename = paste0(
#             figdir, '/', filename_base, 'UMAP_metadata_',
#             sample
#         ),
#         width = 6,
#         height = 5,
#         res = 300,
#         png_only = TRUE
#     )

# }



# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sample_ID"

# p = if (group_var == split_var) {
#     p = CellDimPlot(seu, group_by = group_var, reduction = "wnn.umap_harmony")
# } else {
#     p = CellDimPlot(seu, split_by = split_var, group_by = group_var, reduction = "wnn.umap_harmony")
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 16,
#     height = 8
# )

# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sex"

# p = if (group_var == split_var) {
#     p = CellDimPlot(seu, group_by = group_var, reduction = "wnn.umap_harmony")
# } else {
#     p = CellDimPlot(seu, split_by = split_var, group_by = group_var, reduction = "wnn.umap_harmony")
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 16,
#     height = 8
# )


# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Condition"

# p = if (group_var == split_var) {
#     p = CellDimPlot(seu, group_by = group_var, reduction = "wnn.umap_harmony")
# } else {
#     p = CellDimPlot(seu, split_by = split_var, group_by = group_var, reduction = "wnn.umap_harmony")
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 16,
#     height = 8
# )


# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Condition"
# split_val = "MelanieRuffner-08182023"

# if (group_var == split_var) {
#     p = CellDimPlot(
#         subset(seu, subset = Dataset == split_val),
#         group_by = group_var,
#         reduction = "wnn.umap_harmony"
#     )
# } else {
#     p = CellDimPlot(
#         subset(seu, subset = Dataset == split_val),
#         split_by = split_var,
#         group_by = group_var,
#         reduction = "wnn.umap_harmony"
#     )
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 8,
#     height = 8
# )

# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sample Status"
# split_val = "HTS_SO243"

# if (group_var == split_var) {
#     p = CellDimPlot(
#         subset(seu, subset = Dataset == split_val),
#         group_by = group_var,
#         reduction = "wnn.umap_harmony"
#     )
# } else {
#     p = CellDimPlot(
#         subset(seu, subset = Dataset == split_val),
#         split_by = split_var,
#         group_by = group_var,
#         reduction = "wnn.umap_harmony"
#     )
# }
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_',
#         if (group_var == split_var) group_var else paste0(group_var, '__', split_var)
#     ),
#     width = 8,
#     height = 8
# )



# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# seu@meta.data[['Condition']] = factor(
#     seu@meta.data[['Condition']],
#     levels = c("Control", 
#     "Remission (no dupilumab)", 
#     "Remission (dupilumab)",
#     "Active (dupilumab)", 
#     "Active (no dupilumab)")
# )
# split_var = "Dataset"
# group_var = "Condition"
# color_val = "Sample_ID"
# split_val = "MelanieRuffner-08182023"


# p = CellDimPlot(
#     subset(seu, subset = Dataset == split_val),
#     split_by = group_var,
#     group_by = color_val,
#     reduction = "wnn.umap_harmony"
# )

# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_',
#         if (color_val == split_var) color_val else paste0(color_val, '__', split_var)
#     ),
#     width = 18,
#     height = 12
# )

# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sample Status"
# color_val = "Sample_ID"
# split_val = "HTS_SO243"


# p = CellDimPlot(
#     subset(seu, subset = Dataset == split_val),
#     split_by = group_var,
#     group_by = color_val,
#     reduction = "wnn.umap_harmony"
# )

# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_',
#         if (color_val == split_var) color_val else paste0(color_val, '__', split_var)
#     ),
#     width = 16,
#     height = 10
# )


# # # --  ---------------------------- ##
# # feature = 'nCount_RNA'
# # p = FeatureStatPlot(
# #     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
# #     plot_type = "violin",
# #     group.by = "Sample_ID",
# #     features = feature,
# #     pt_size = 0.001
# # ) + NoLegend()
# # plot_pdf_and_png(
# #     plot = p, 
# #     filename = paste0(
# #         figdir, '/', filename_base, 'VlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
# #     ),
# #     width = 16,
# #     height = 6
# # )

# # # --  ---------------------------- ##
# # feature = 'nCount_RNA'
# # p = FeatureStatPlot(
# #     subset(seu, subset = Dataset == "HTS_SO243"),
# #     plot_type = "violin",
# #     group.by = "Sample_ID",
# #     features = feature,
# #     pt.size = 0.001
# # ) + NoLegend()
# # plot_pdf_and_png(
# #     plot = p, 
# #     filename = paste0(
# #         figdir, '/', filename_base, 'VlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
# #     ),
# #     width = 16,
# #     height = 6
# # )







# # -- UMAP: Sample_ID x Dataset ---------------------------- ##
# split_var = "Dataset"
# group_var = "Sample_ID"
# split_val = "HTS_SO243"


# p = CellDimPlot(
#     subset(seu, subset = Dataset == split_val),
#     group_by = group_var,
#     reduction = "wnn.umap_harmony"
# )

# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_', group_var
#     ),
#     width = 12,
#     height = 10
# )



# split_var = "Dataset"
# group_var = "Sample_ID"
# split_val = "MelanieRuffner-08182023"


# p = CellDimPlot(
#     subset(seu, subset = Dataset == split_val),
#     group_by = group_var,
#     reduction = "wnn.umap_harmony"
# )

# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterUMAP_metadata_', split_val, '_', group_var
#     ),
#     width = 12,
#     height = 10
# )























# # --  ---------------------------- ##
# feature = 'nCount_RNA'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'nFeature_RNA'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'nCount_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'nFeature_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'percent_mt'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'IsoCount'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )





# ###################################


# # --  ---------------------------- ##
# feature = 'nCount_RNA'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )

# # --  ---------------------------- ##
# feature = 'nFeature_RNA'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )

# # --  ---------------------------- ##
# feature = 'nCount_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )

# # --  ---------------------------- ##
# feature = 'nFeature_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )

# # --  ---------------------------- ##
# feature = 'percent_mt'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )

# # --  ---------------------------- ##
# feature = 'IsoCount'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "ridge",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterRidgePlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 10,
#     height = 12
# )


# ###################################


















# # --  ---------------------------- ##
# feature = 'nCount_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'nFeature_ADT'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'percent_mt'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )

# # --  ---------------------------- ##
# feature = 'IsoCount'
# p = FeatureStatPlot(
#     seu,
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', '_', feature, '_', "Sample_ID"
#     ),
#     width = 20,
#     height = 8
# )


# # seurat = FindClusters(
# # 	seu,
# # 	algorithm = 3,
# # 	graph.name = 'wsnn_harmony',
# # 	resolution = c(seq(0, 1, 0.2))
# # )













































# # --  ---------------------------- ##
# feature = 'nCount_RNA'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )


# # --  ---------------------------- ##
# feature = 'nFeature_RNA'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )

# # --  ---------------------------- ##
# feature = 'nCount_ADT'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )


# # --  ---------------------------- ##
# feature = 'nFeature_ADT'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )


# # --  ---------------------------- ##
# feature = 'IsoCount'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )

# # --  ---------------------------- ##
# feature = 'percent_mt'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "MelanieRuffner-08182023"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt.size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "MelanieRuffner-08182023", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )









# # --  ---------------------------- ##
# feature = 'nCount_RNA'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )

# # --  ---------------------------- ##
# feature = 'nFeature_RNA'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )
# # --  ---------------------------- ##
# feature = 'nCount_ADT'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )

# # --  ---------------------------- ##
# feature = 'nFeature_ADT'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )
# # --  ---------------------------- ##
# feature = 'percent_mt'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )
# # --  ---------------------------- ##
# feature = 'IsoCount'
# p = FeatureStatPlot(
#     subset(seu, subset = Dataset == "HTS_SO243"),
#     plot_type = "violin",
#     group_by = "Sample_ID",
#     features = feature,
#     pt_size = 0.001,
#     add_point = TRUE
# ) + NoLegend()
# plot_pdf_and_png(
#     plot = p, 
#     filename = paste0(
#         figdir, '/', filename_base, 'scplotterVlnPlot_metadata_', "HTS_SO243", '_', feature, '_', "Sample_ID"
#     ),
#     width = 16,
#     height = 6
# )



# table(seu@meta.data$Sample_ID) %>% as.data.frame()
# write.table(
#     table(seu@meta.data$Sample_ID) %>% as.data.frame(),
#     file = 'cell_counts.tsv',
#     sep = '\t',
#     row.names = FALSE,
#     quote = FALSE
# )
