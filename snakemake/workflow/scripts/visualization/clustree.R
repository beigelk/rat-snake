
# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_obj      = snakemake@input[['seurat_obj']]

# params
project           = snakemake@params[['project']]
clustree_json     = snakemake@params[['clustree_json']]
figbase           = snakemake@params[['figbase']]
tblbase           = snakemake@params[['tblbase']]
r_library         = snakemake@params[['r_library']]

# output
seurat_obj_clustered = snakemake@output[['seurat_obj_clustered']]

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

# utility
library(qs2)
library(tidyverse)
library(stringdist)
library(jsonlite)


# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
# options(future.globals.maxSize = 200000 * 1024^2)


# -- SOURCE SCRIPTS ------------------------- #
sapply(c('snakemake/workflow/scripts/utils/plot.R'), source, print.eval = FALSE)


# -- FUNCTIONS ------------------------------ #


# -- MAKE OUPUT DIRS ------------------------ #
figdir = dirname(figbase)
tbldir = dirname(tblbase)
objdir = dirname(seurat_obj_anno)
sapply(
    c(figdir, tbldir, objdir),
    dir.create, showWarnings = FALSE, recursive = TRUE
)


# == MAIN =================================== #
# -- LOAD OBJECT ---------------------------- #
seu = qs2::qs_read(
    file = seurat_obj
)

# -- CLUSTREE ------------------------------- #
dir.create(paste0(figdir, 'clustree'), showWarnings = FALSE, recursive = TRUE)
clustree = clustree(
    seu_lympho,
    prefix = "wsnn_harmony_res.",
    ) +
    ggtitle(project) +
    theme(
        legend.spacing = unit(0.1, "cm"), 
        legend.key.size = unit(0.7, 'lines')
    )

plot_pdf_and_png(
    plot = clustree,
    filename = paste0(figbase),
    width = 12,
    height = 8
)


# -- UMAPS OF CLUSTERS ---------------------- #
# Plot dimplots of all of the resolutions
# --------------------------------------------------------------------------
for (col in colnames(seu_lympho@meta.data)[grepl('wsnn_harmony_res*', colnames(seu_lympho@meta.data))]){

    dir.create(paste0(figdir, 'clusters/'), showWarnings = FALSE, recursive = TRUE)
    plot_pdf_and_png(
        plot = DimPlot(seu_lympho, group.by = col, label = TRUE, reduction = 'wnn.umap_harmony'),
        filename = paste0(figdir, 'clusters/', project, '__DimPlot__', col),
        height = 5, 
        width = 5
    )
}
# --------------------------------------------------------------------------

# Get lists of markers
# --------------------------------------------------------------------------
marker_list = list(
    adt_lymphocyte_markers = read.table(marker_file)$V1,
    markers_additional = read.table(marker_file_2, sep = '\t')$V1
)

actual_adts = rownames(seu_lympho@assays$ADT$counts)

# String matching to make sure we have the correct
normalize_marker <- function(name) {
  name <- tolower(name)
  name <- gsub("-", "", name)
  name <- gsub("α", "alpha", name)
  name <- gsub("β", "beta", name)
  name <- gsub("δ", "delta", name)
  name <- gsub("\\s+", "", name)
  name <- gsub("\\(|\\)", "", name)
  return(name)
}

# Normalize both
norm_markers <- sapply(marker_list[[2]], normalize_marker)
norm_columns <- sapply(actual_adts, normalize_marker)

# --- 2. Matching function ---
match_markers <- function(query_markers, dataset_columns, maxDist = 0.4) {

  norm_query <- sapply(query_markers, normalize_marker, USE.NAMES = FALSE)
  norm_cols  <- sapply(dataset_columns, normalize_marker, USE.NAMES = FALSE)
  results <- character(length(query_markers))

  for (i in seq_along(norm_query)) {
    q <- norm_query[i]
    q_raw <- query_markers[i]

    # Extract CD number (if present)
    cd_num <- sub(".*cd([0-9]+).*", "\\1", q)
    has_cd_num <- grepl("cd[0-9]+", q)

    # --- Step 1: exact match ---
    exact_idx <- which(norm_cols == q)
    if (length(exact_idx) == 1) {
      results[i] <- dataset_columns[exact_idx]
      next
    }

    # --- Step 2: word boundary match (e.g., cd22 should not match cd223) ---
    word_idx <- grep(paste0("\\b", q, "\\b"), norm_cols, perl = TRUE)
    if (length(word_idx) == 1) {
      results[i] <- dataset_columns[word_idx]
      next
    }

    # --- Step 3: partial substring (safe numeric boundaries) ---
    part_idx <- grep(paste0("(^|[^0-9])", q, "([^0-9]|$)"), norm_cols, perl = TRUE)
    if (length(part_idx) == 1) {
      results[i] <- dataset_columns[part_idx]
      next
    }

    # --- Step 4: restrict fuzzy candidates ---
    # if query is a CD marker, only compare to same CD number
    candidate_idx <- if (has_cd_num) {
      grep(paste0("cd", cd_num, "(?![0-9])"), norm_cols, perl = TRUE)
    } else seq_along(norm_cols)

    # --- Step 5: fuzzy match only within filtered candidates ---
    if (length(candidate_idx) > 0) {
      dists <- stringdist(q, norm_cols[candidate_idx], method = "jw")
      min_dist <- min(dists)
      if (min_dist < maxDist) {
        best_idx <- candidate_idx[which.min(dists)]
        results[i] <- dataset_columns[best_idx]
        next
      }
    }

    # --- Step 6: no match found ---
    results[i] <- NA
  }

  names(results) <- query_markers
  results
}

matches = match_markers(marker_list[[2]], actual_adts)

# Update marker list with the correct ADT names
marker_list[[2]] = unique(unname(matches))
# --------------------------------------------------------------------------

# For all resolutions, make plots of markers
# --------------------------------------------------------------------------
DefaultAssay(seu_lympho) = 'ADT'

# DOTPLOTS
for (item in 1:length(marker_list)) {

    # Make a subdir for the marker list
    marker_setname = names(marker_list)[item]
    markers = marker_list[[item]]

    figdir_dotplot = paste0(figdir, 'adt/', marker_setname, '/dotplots/')
    for (subdir in c(figdir_dotplot)) {
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    }

    res_list = colnames(seu_lympho@meta.data)[grep('wsnn_harmony', colnames(seu_lympho@meta.data))]

    # For each resolution, 
    for (res in res_list) {

        dotplot = DotPlot(
            seu_lympho,
            features = markers,
            group.by = res
        ) +
        labs(title = paste0(project, ' ', res, ': ', marker_setname)) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        )

        plot_pdf_and_png(
            plot = dotplot,
            filename = paste0(figdir_dotplot, project, '__DotPlot__', marker_setname, '__', res),
            width = 12,
            height = 8
        ) 
    }
}

# CLUSTREE PLOTS
for (item in 1:length(marker_list)) {

    # Make a subdir for the marker list
    marker_setname = names(marker_list)[item]
    markers = marker_list[[item]]

    figdir_clustree = paste0(figdir, 'adt/', marker_setname, '/clustree/')
    for (subdir in c(figdir_clustree)) {
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    }

    res_list = colnames(seu_lympho@meta.data)[grep('wsnn_harmony', colnames(seu_lympho@meta.data))]

    for (marker in markers) {

        print(marker)

        # CLUSTREE
        clustree = clustree(
            seu_lympho,
            prefix = "wsnn_harmony_res.",
            node_colour = marker,
            node_colour_aggr = "median",
            node_text_colour = "white"
            ) + 
            ggtitle(paste0('EoE v4.7 Lymphocytes: ', marker))

        plot_pdf_and_png(
            plot = clustree,
            filename = paste0(figdir_clustree, project, '__Clustree__', marker_setname, '__', marker),
            width = 12,
            height = 8
        )
    }
}

# FEATURE PLOTS
for (item in 1:length(marker_list)) {

    # Make a subdir for the marker list
    marker_setname = names(marker_list)[item]
    markers = marker_list[[item]]

    figdir_featplot = paste0(figdir, 'adt/', marker_setname, '/featplots/')
    for (subdir in c(figdir_featplot)) {
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    }

    for (marker in markers) {

        print(marker)

        # UMAP
        umap = FeaturePlot(seu_lympho, features = marker, reduction = 'wnn.umap_harmony') +
            ggtitle(paste0('EoE v4.7 Lymphocytes: ', marker))
    
        plot_pdf_and_png(
            plot = umap,
            filename = paste0(figdir_featplot, project, '__FeaturePlot__', marker_setname, '__', marker),
            width = 5,
            height = 5
        )
    }

}
# -----------------------------------------------------------------------------------------------

# RIDGE PLOTS
# -----------------------------------------------------------------------------------------------
for (item in 1:length(marker_list)) {

    # Make a subdir for the marker list
    marker_setname = names(marker_list)[item]
    markers = marker_list[[item]]

    figdir_ridge = paste0(figdir, 'adt/', marker_setname, '/ridgeplots/')
    for (subdir in c(figdir_ridge)) {
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    }

    res_list = colnames(seu_lympho@meta.data)[grep('wsnn_harmony', colnames(seu_lympho@meta.data))]

    # For each resolution, 
    for (res in res_list) {

        h = 7
        h = h +
            (seu_lympho@meta.data[, res] %>%
                levels() %>%
                length()
            )/5

        # ridge = RidgePlot(seu_lympho, group.by = res, features = markers, stack = TRUE) +
        #     ggtitle(paste0('EoE v4.7 Lymphocytes: ', marker_setname)) +
        #     NoLegend()
        ridge = RidgePlot(
            seu_lympho, features = markers, group.by = res, stack = TRUE, fill.by = 'ident', same.y.lims = TRUE
        ) +
        theme(
            plot.title = element_text(size = 12),
            strip.text.y.right = element_text(size = 10, hjust = 0)
        ) +
        labs(
                title = paste0(project, ' ', res, ': ', marker_setname)
        ) +
        NoLegend()
    
        plot_pdf_and_png(
            plot = ridge,
            filename = paste0(figdir_ridge, project, '__RidgePlot__', res, '__', marker_setname),
            width = 20,
            height = h
        )
    }
}
# -----------------------------------------------------------------------------------------------









# FindAllMarkers at Resolution 1
# -----------------------------------------------------------------------------------------------
source('Rscripts/utils.R')

seu_lympho = qs2::qs_read(
    file = seurat_outfile
)
# -----------------------------------------------------------------------------------------------

# RNA
# -----------------------------------------------------------------------------------------------
outdir = 'snakemake/EoE_Dupilumab_V4/output/v4.7.0.0/annotation/lymphocytes/tables/FindAllMarkers/RNA/'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# colnames(seu_lympho@meta.data)
# DimPlot(seu_lympho, group.by = 'wsnn_harmony_res.1')
res = 'wsnn_harmony_res.1'

DefaultAssay(seu_lympho) = 'RNA'
Idents(seu_lympho) = res
all_markers = FindAllMarkers(
    seu_lympho,
    group.by = res
)

# Save full table to CSV
write.table(
    all_markers,
    file = paste0(outdir, project, '__RNA__', res, '__', 'FindAllMarkers__all_clusters.tsv'),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Get unique cluster IDs
clusters = unique(all_markers$cluster)

# Loop through clusters
for (clust in clusters) {
    cluster_markers = all_markers %>% filter(cluster == clust)
    filename = paste0(outdir, project, '__RNA__', res, '__', 'FindAllMarkers__cluster_', clust, '.tsv')
    write.table(
        cluster_markers,
        file = filename,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

all_markers %>%
    group_by(cluster) %>%
    dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 0.5
    ) %>%
    slice_head(n = 5) %>%
    ungroup() -> top5

all_markers %>%
    group_by(cluster) %>%
    dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 0.5
    ) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

seurat_subset = seu_lympho
seurat_subset = ScaleData(seurat_subset, features = top10$gene)

plot_pdf_and_png(
    plot = DoHeatmap(seurat_subset, features = top5$gene) + NoLegend() + labs(title = paste0(project, ' ', res, ': cluster ', clust), subtitle = 'RNA data'),

    file = paste0(outdir, project, '__RNA__', res, '__', 'FindAllMarkers__Heatmap'),
    width = 18,
    height = 12
)

plot_pdf_and_png(
    plot = DotPlot(seurat_subset, features = unique(top5$gene)) +
                theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(title = paste0(project, ' ', res, ': cluster ', clust), subtitle = 'RNA data'),

    file = paste0(outdir, project, '__RNA__', res, '__', 'FindAllMarkers__DotPlot'),
    width = 20,
    height = 8
)


for (clust in clusters) {

    plot_pdf_and_png(
        plot = VlnPlot(
            seurat_subset,
            features = unique(top10 %>% filter(cluster == clust) %>% pull(gene)),
            stack = TRUE,
            flip = TRUE
        ) +
            theme(
                strip.text.y.right = element_text(hjust = 0)
        ) +
            NoLegend() +
            labs(title = paste0(project, ' ', res, ': cluster ', clust),
            subtitle = 'RNA data'
        ),
        file = paste0(outdir, project, '__RNA__', res, '__', 'FindAllMarkers__VlnPlot__cluster_', clust),
        width = 10,
        height = 8
    )
}
# -----------------------------------------------------------------------------------------------


# ADT
# -----------------------------------------------------------------------------------------------
outdir = 'snakemake/EoE_Dupilumab_V4/output/v4.7.0.0/annotation/lymphocytes/tables/FindAllMarkers/ADT/'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# colnames(seu_lympho@meta.data)
# DimPlot(seu_lympho, group.by = 'wsnn_harmony_res.1')
res = 'wsnn_harmony_res.1'

DefaultAssay(seu_lympho) = 'ADT'
Idents(seu_lympho) = res
all_markers = FindAllMarkers(
    seu_lympho,
    group.by = res
)

# Save full table to CSV
write.table(
    all_markers,
    file = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__all_clusters.tsv'),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Get unique cluster IDs
clusters = unique(all_markers$cluster)

# Loop through clusters
for (clust in clusters) {
    cluster_markers = all_markers %>% filter(cluster == clust)
    filename = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__cluster_', clust, '.tsv')
    write.table(
        cluster_markers,
        file = filename,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

all_markers %>%
    group_by(cluster) %>%
    dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 0.5
    ) %>%
    slice_head(n = 5) %>%
    ungroup() -> top5

all_markers %>%
    group_by(cluster) %>%
    dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 0.5
    ) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

seurat_subset = seu_lympho
seurat_subset = ScaleData(seurat_subset, features = top10$gene)

plot_pdf_and_png(
    plot = DoHeatmap(seurat_subset, features = top5$gene) + NoLegend() + labs(title = paste0(project, ' ', res, ': cluster ', clust), subtitle = 'ADT data'),

    file = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__Heatmap'),
    width = 18,
    height = 12
)

plot_pdf_and_png(
    plot = DotPlot(seurat_subset, features = unique(top5$gene)) +
                theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                labs(
                    title = paste0(project, ' ', res, ': cluster ', clust),
                    subtitle = 'ADT data'
            ),

    file = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__DotPlot'),
    width = 20,
    height = 8
)

for (clust in clusters) {

    if (length(unique(top10 %>% filter(cluster == clust) %>% pull(gene))) > 1) {

    #     plot_pdf_and_png(
    #         plot = VlnPlot(
    #             seurat_subset,
    #             features = unique(top10 %>% filter(cluster == clust) %>% pull(gene)),
    #             stack = TRUE,
    #             flip = TRUE
    #         ) +
    #             theme(
    #                 strip.text.y.right = element_text(hjust = 0)
    #         ) +
    #             NoLegend() +
    #             labs(
    #                 title = paste0(project, ' ', res, ': cluster ', clust),
    #                 subtitle = 'ADT data'
    #         ),

    #     file = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__VlnPlot__cluster_', clust),
    #     width = 10,
    #     height = 8
    # )

    } else if (length(unique(top10 %>% filter(cluster == clust) %>% pull(gene))) == 1) {
        plot_pdf_and_png(
            plot = VlnPlot(
                seurat_subset,
                features = unique(top10 %>% filter(cluster == clust) %>% pull(gene)),
                pt.size = 0,
                cols = rep('#F8766D', length(unique(seurat_subset@meta.data[[res]])))
            ) +
            #     theme(
            #         strip.text.y.right = element_text(hjust = 0)
            # ) +
                NoLegend() +
                labs(
                    title = paste0(project, ' ', res, ': cluster ', clust),
                    subtitle = 'ADT data'
            ),

        file = paste0(outdir, project, '__ADT__', res, '__', 'FindAllMarkers__VlnPlot__cluster_', clust),
        width = 10,
        height = 5
    )
    }


}
# -----------------------------------------------------------------------------------------------

