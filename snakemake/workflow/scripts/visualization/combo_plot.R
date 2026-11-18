
# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj = snakemake@input[['seurat_obj']]

# params
assay       = snakemake@params[['assay']]
marker_dirs = snakemake@params[['marker_dirs']]
res_prefix  = snakemake@params[['res_prefix']]
outbase     = snakemake@params[['outbase']]
r_library   = snakemake@params[['r_library']]

# output
done = snakemake@output[['done']]

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
library(qs2)
library(tidyverse)

# plotting
library(ggplot2)
library(clustree)
library(scplotter)
library(patchwork)
library(plotthis)

# analysis
library(Seurat)


# -- OUTPUT DIR --------------------------------------- #
dir.create('V5', showWarnings = FALSE, recursive = TRUE)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
# plan("multicore", workers = 16)
# options(future.globals.maxSize = 200000 * 1024^2)


# -- SOURCE SCRIPTS ------------------------- #
sapply(
    paste0('snakemake/workflow/scripts/utils/', c('file_io.R', 'plot.R')),
    source,
    print.eval = FALSE
)


# -- FUNCTIONS ------------------------------ #


# -- LOAD OBJECT ---------------------------- #
seu = qs2::qs_read(
    seurat_obj
)


# -- GET MARKER LISTS ----------------------- #
marker_lists = sapply(
    marker_dirs,
    list.files
)


# -- BASIC FEATURE PLOTS -------------------- #
for (ml in 1:length(marker_lists)) {
    
    group   = basename(names(marker_lists)[ml])
    ml_files = file.path(names(marker_lists)[ml], marker_lists[[ml]])

    for (f in ml_files) {    
        
        ml_name = paste0(group, '__', tools::file_path_sans_ext(basename(f)))

        # Outfile name
        outfile = paste0(outbase, ml_name, "__FeatureDim")

        # Load marker list
        markers = read.table(
            file = f,
            col.names = ml_name,
            check.names = FALSE
        ) %>%
        pull(ml_name)

        # plotthis::show_palettes(type = "continuous")
        # SunsetDark
        # ocean.dense
        pal = "ocean.dense"
        p = FeatureStatPlot(
            seu,
            assay = assay,
            title = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
            plot_type = "dim",
            features = markers,
            reduction = "wnn.umap_harmony",
            bg_cutoff = -Inf,
            # pt_size = 1,
            palette = pal,
            palreverse = FALSE,
            facet_ncol = 3,
            raster = FALSE
        ) +
        theme(
            plot.title = element_text(size = 24),
            strip.text.x = element_text(size = 18)
        )

        png(
            filename = paste0(outfile, ".png"),
            width = 15,
            height = 18 * length(markers) / 10,
            units = "in",
            res = 600
        )
        print(p)
        dev.off()
    }
}

# -- STACK VIOLIN PLOTS -------------------- #
for (ml in 1:length(marker_lists)) {
    
    group   = basename(names(marker_lists)[ml])
    ml_files = file.path(names(marker_lists)[ml], marker_lists[[ml]])

    for (f in ml_files) {    
        
        ml_name = paste0(group, '__', tools::file_path_sans_ext(basename(f)))

        res_list = colnames(seu@meta.data)[grep(res_prefix, colnames(seu@meta.data))]

        # Load marker list
        markers = read.table(
            file = f,
            col.names = ml_name,
            check.names = FALSE
        ) %>%
        pull(ml_name)

        for (res in res_list) {
            
            # Outfile name
            res_val = gsub(res_prefix, '', res)
            outfile = paste0(outbase, ml_name, "__StackViolin__res_", res)

            n_clust = length(unique(seu@meta.data[[res]]))

            # plotthis::show_palettes(type = "continuous")
            # SunsetDark
            # ocean.dense
            # pal = "ocean.solar"
            pal = "jcolors-rainbow"
            p = FeatureStatPlot(
                seu,
                assay = assay,
                title = paste0("Resolution ", res_val, " (", res, ")"),
                subtitle = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
                plot_type = "violin",
                features = markers,
                ident = res,
                add_bg = TRUE,
                add_box = TRUE,
                box_ptsize = 1,
                stack = TRUE,
                palette = pal,
                palreverse = FALSE,
                raster = FALSE
            )

            p = p +
                theme(
                    plot.title    = element_text(size = 14),
                    strip.text.x  = element_text(size = 14),
                    plot.subtitle = element_text(size = 12)
                ) + 
                NoLegend()

            png(
                filename = paste0(outfile, ".png"),
                width = 3.5 + (n_clust / 2),
                height = 14 * length(markers) / 10,
                units = "in",
                res = 600
            )
            print(p)
            dev.off()
        }
    }
}

# res_list = c(0.2, 0.4, 0.6, 0.8)
# -- COMBO PLOTS ------------------------------ #
for (ml in 1:length(marker_lists)) {
    
    group   = basename(names(marker_lists)[ml])
    ml_files = file.path(names(marker_lists)[ml], marker_lists[[ml]])

    for (f in ml_files) {    
        
        ml_name = paste0(group, '__', tools::file_path_sans_ext(basename(f)))
        res_prefix = 'cluster__epithelial__wsnn_harmony_res.'
        res_list = unique(trimws(colnames(seu@meta.data)[grep(res_prefix, colnames(seu@meta.data))]))

        # Load marker list
        markers = read.table(
            file = f,
            col.names = ml_name,
            check.names = FALSE
        ) %>%
        pull(ml_name)

        feat_ncol = if (length(markers) < 11) 5 else 12
        row1_margin = if (length(markers) < 11) 60 else -10

        pal = "ocean.dense"
        p3 = FeatureStatPlot(
            seu,
            assay = assay,
            # title = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
            plot_type = "dim",
            features = markers,
            reduction = "wnn.umap_harmony",
            bg_cutoff = -Inf,
            # pt_size = 1,
            palette = pal,
            palreverse = TRUE,
            facet_ncol = feat_ncol,
            # raster = FALSE
            raster = FALSE
        ) + 
        guides(
            colour = guide_colorbar(title = "Expression")
        ) +
        theme(
            plot.margin = margin(t = row1_margin, b = row1_margin),
            legend.key.height = unit(0.8, "cm"),
            legend.key.width  = unit(0.3, "cm"),
            legend.margin = margin(0, 0, 0, -3)
        )

        for (res in res_list) {
            
            # Outfile name
            res_val = gsub(res_prefix, '', res)
            outfile = paste0(outbase, '/', ml_name, "__ComboPlot__res_", res)

            n_clust = length(unique(seu@meta.data[[res]]))
            res_col = paste0('Cluster')

            seu@meta.data[[res_col]] = paste0(seu@meta.data[[res]], ' ')
            seu@meta.data[[res_col]] = factor(
                seu@meta.data[[res_col]],
                levels = paste0(levels(seu@meta.data[[res]]), ' ')
            )

            row2_margin = if (length(markers) < 11) -40 else -80
            flip_vln = if (length(markers) < 11) FALSE else TRUE

            # -- VIOLIN
            # plotthis::show_palettes(type = "continuous")
            # SunsetDark
            # ocean.dense
            # pal = "ocean.solar"
            row2_t_margin = -30
            row2_margin = -40
            pal = "jcolors-rainbow"
            p1 = FeatureStatPlot(
                seu,
                assay = assay,
                # title = paste0("Resolution ", res_val, " (", res, ")"),
                # subtitle = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
                plot_type = "violin",
                features = markers,
                ident = res_col,
                add_bg = TRUE,
                add_box = TRUE,
                box_ptsize = 1,
                stack = TRUE,
                flip = flip_vln,
                palette = pal,
                palreverse = FALSE,
                raster = FALSE
            ) +
            theme(
                plot.title    = element_text(size = 14),
                strip.text.x  = element_text(size = 14),
                plot.subtitle = element_text(size = 12),
                plot.margin = margin(t = row2_t_margin, b = row2_margin)
            ) + 
            NoLegend()

            # -- UMAP
            p2 = CellDimPlot(
                seu,
                # title = paste0("Resolution ", res_val, " (", res, ")"),
                group_by = res_col,
                label = TRUE,
                label_insitu = TRUE,
                label_size = 5,
                label_repel = TRUE,
                label_repulsion = 10,
                label_bg_r = 0.001,
                label_bg = "white",
                label_pt_size = 0.1,
                label_pt_color = "#00000000", 
                label_segment_color = "#00000000", 
                label_fg = "black",
                palette = pal,
                palreverse = FALSE,
                raster = FALSE
            ) +
            theme(
                plot.title    = element_text(size = 14),
                plot.margin   = margin(t = row2_t_margin, b = row2_t_margin, l = 20)
            )
            
            # -- WRITE PLOT
            h = if (length(markers) < 11) 14 else 3 + (length(markers) / 2)
            w = if (length(markers) < 11) 16 else 6 + (length(markers))
            ratios = if (length(markers) < 11) c(6, 0.01, 5) else c(7 + length(markers) / 2, -2, 6 + length(markers) / 2)
            
            png(
                filename = basename(paste0(outfile, ".png")),
                width = w,
                height = h,
                # width = 3.5 + (n_clust / 2),
                # height = 14 * length(markers) / 10,
                units = "in",
                res = 600
            )
            # print((p1 | p2) / p3)
            print(
                (p1 | p2) / plot_spacer() / p3 +
                plot_layout(heights = ratios) +
                plot_annotation(
                    title = paste0("Resolution ", res_val, " (", res, ")"),
                    subtitle = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
                    theme = theme(
                        plot.title = element_text(
                            size = 24,
                            face = "bold",
                            hjust = 0.5
                        ),
                        plot.subtitle = element_text(
                            size = 20,
                            hjust = 0.5
                        )
                    )
                ) 
                 #+
                # plot_layout(heights = c(1, 2))
            )
            dev.off()
        }
    }
}


# # -- HEATMAP PLOTS --------------------------- #
# for (ml in 1:length(marker_lists)) {
    
#     group   = basename(names(marker_lists)[ml])
#     ml_files = file.path(names(marker_lists)[ml], marker_lists[[ml]])

#     rows_data = data.frame(
#         Features     = character(),
#         `Cell types` = character()
#     )

#     for (f in ml_files) {
        
#         ml_name = paste0(group, '__', tools::file_path_sans_ext(basename(f)))

#         res_list = colnames(seu@meta.data)[grep(res_prefix, colnames(seu@meta.data))]

#         # Load marker list
#         markers = read.table(
#             file = f,
#             col.names = ml_name,
#             check.names = FALSE
#         ) %>%
#         pull(ml_name)

#         rows_data = bind_rows(
#             rows_data,
#             data.frame(
#                 Features     = markers,
#                 `Cell types` = rep(ml_name, length(markers))
#             )
#         )
#     }

#     for (res in res_list) {
        
#         # Outfile name
#         res_val = gsub(res_prefix, '', res)
#         outfile = paste0(outbase, ml_name, "__StackViolin__res_", res)

#         n_clust = length(unique(seu@meta.data[[res]]))

#         # plotthis::show_palettes(type = "continuous")
#         # SunsetDark
#         # ocean.dense
#         # pal = "ocean.solar"
#         pal = "jcolors-rainbow"
#         p = FeatureStatPlot(
#             seu,
#             assay = assay,
#             title = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
#             subtitle = paste0("Resolution ", res_val, " (", res, ")"),
#             plot_type = "heatmap",
#             features = markers,
#             ident = res,
#             cell_type = "dot",
#             name = "Expression Level",
#             dot_size = function(x) sum(x > 0) / length(x),
#             dot_size_name = "Percent Expressed",
#             add_bg = TRUE,
#             rows_data = rows_data,
#             # show_row_names = TRUE,
#             # rows_split_by = "group",
#             # cluster_rows = FALSE,
#             column_annotation = c("Condition", "nFeature_RNA", "nFeature_ADT"),
#             column_annotation_type = list(Condition = "pie", nFeature_RNA = "violin", nFeature_ADT = "violin"),
#             column_annotation_params = list(nFeature_RNA = list(show_legend = FALSE), nFeature_ADT = list(show_legend = FALSE)),
#             # row_annotation_side = "right",
#             # row_annotation = c("TF", "CSPA"),
#             # row_annotation_type = list(TF = "simple", CSPA = "simple"),
#             flip = TRUE,
#             raster = FALSE
#         )

#         p = p +
#             theme(
#                 plot.title    = element_text(size = 14),
#                 strip.text.x  = element_text(size = 12),
#                 plot.subtitle = element_text(size = 14)
#             ) + 
#             NoLegend()

#         png(
#             filename = paste0(outfile, ".png"),
#             width = 3.5 + (n_clust / 2),
#             height = 14 * length(markers) / 10,
#             units = "in",
#             res = 600
#         )
#         print(p)
#         dev.off()
#     }
# }



# # -- UMAP
# p2 = CellDimPlot(
#     seu,
#     # title = paste0("Resolution ", res_val, " (", res, ")"),
#     group_by = res_col,
#     split_by = res_col,
#     label = TRUE,
#     label_insitu = TRUE,
#     label_size = 5,
#     label_repel = TRUE,
#     label_repulsion = 10,
#     label_bg_r = 0.001,
#     label_bg = "white",
#     label_pt_size = 0.1,
#     label_pt_color = "#00000000", 
#     label_segment_color = "#00000000", 
#     label_fg = "black",
#     palette = pal,
#     palreverse = FALSE,
#     raster = FALSE
# ) +
# theme(
#     plot.title    = element_text(size = 14),
#     plot.margin   = margin(t = row2_t_margin, b = row2_t_margin, l = 20)
# )






# FeatureStatPlot(
#             seu,
#             assay = assay,
#             # title = toupper(gsub('_', ' ', gsub('__', ': ', ml_name))),
#             split_by = res_col,
#             plot_type = "dim",
#             features = markers,
#             reduction = "wnn.umap_harmony",
#             bg_cutoff = -Inf,
#             # pt_size = 1,
#             palette = pal,
#             palreverse = TRUE,
#             facet_ncol = feat_ncol,
#             # raster = FALSE
#             raster = TRUE
#         ) + 
#         guides(
#             colour = guide_colorbar(title = "Expression")
#         ) +
#         theme(
#             plot.margin = margin(t = row1_margin, b = row1_margin),
#             legend.key.height = unit(0.8, "cm"),
#             legend.key.width  = unit(0.3, "cm"),
#             legend.margin = margin(0, 0, 0, -3)
#         )
