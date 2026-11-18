###############################################################################
# Rule:     dsb
# Purpose:  Normalise and denoise one sample's ADT counts with dsb
#
# Description:
#   Unchanged analysis: empty drops are background barcodes with
#   log10 ADT counts in (1.5, 3) and log10 RNA counts < 2.5; cells within
#   3 MADs of the median RNA and ADT library sizes with mitochondrial
#   proportion < 0.14 are kept; DSBNormalizeProtein with isotype controls;
#   then CITE clusters, UMAP and a heatmap of median dsb values per cluster.
#   Every threshold is now an option under preprocess.dsb (same defaults).
#
# Inputs:   snakemake@params$raw_matrix, $filtered_matrix (Cell Ranger multi)
# Outputs:  adt_dsb (qs2 list: seurat, cells.dsb.norm, cellmd), adt_scat_png,
#           adt_heatmap_png, summary (TSV: sample, dsb_cells, background_drops)
###############################################################################

source(file.path(Sys.getenv("SNAKE_CORE"), "R", "snake_init.R"))
smk <- snake_init(snakemake)

suppressPackageStartupMessages({
  library(dsb)
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(ggplot2)
})


#' Normalise and denoise ADT counts with dsb
#'
#' @param raw.path,filt.path Character. Cell Ranger raw and filtered matrix folders.
#' @param opts List. preprocess.dsb options (thresholds and isotype pattern).
#' @param mt.pattern Character. Regex for mitochondrial genes.
#' @param scat.png,heatmap.png Character. Output figures.
#'
#' @return List: `seurat`, `cells.dsb.norm` (DSBNormalizeProtein output with
#'   stats), `cellmd` (droplet QC of kept cells), `n_background`.
run_dsb <- function(raw.path, filt.path, opts, mt.pattern, scat.png, heatmap.png) {  
    # read raw data using the Seurat function "Read10X_h5"
    raw = Seurat::Read10X(data.dir = raw.path)
    cells = Seurat::Read10X(data.dir = filt.path)

    # define cell-containing barcodes and separate cells and empty drops
    stained_cells = colnames(cells$`Gene Expression`)
    background = setdiff(colnames(raw$`Gene Expression`), stained_cells)

    # split the data into separate matrices for RNA and ADT
    prot = raw$`Antibody Capture`
    rna = raw$`Gene Expression`

    # create metadata of droplet QC stats used in standard scRNAseq processing
    mtgene = grep(pattern = mt.pattern, rownames(rna), value = TRUE) # used below

    md = data.frame(
        rna.size = log10(Matrix::colSums(rna)), 
        prot.size = log10(Matrix::colSums(prot)), 
        n.gene = Matrix::colSums(rna > 0), 
        mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
    )

    # add indicator for barcodes Cell Ranger called as cells
    md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

    # remove barcodes with no evidence of capture in the experiment
    md = md[md$rna.size > 0 & md$prot.size > 0, ]

    png(file = scat.png, height = 6, width = 10, res = 300, units = 'in')
    print(
        ggplot(md, aes(x = log10(n.gene), y = prot.size)) +
        theme_bw() + 
        geom_bin2d(bins = 300) + 
        scale_fill_viridis_c(option = "C") + 
        facet_wrap(~drop.class)
    )
    dev.off()

    background_drops = rownames(
        md[ md$prot.size > opts$background_prot_min &
            md$prot.size < opts$background_prot_max &
            md$rna.size < opts$background_rna_max, ]
    ) 
    background.adt.mtx = as.matrix(prot[ , background_drops])

    # calculate statistical thresholds for droplet filtering.
    cellmd = md[md$drop.class == 'cell', ]

    # filter drops with + / - 3 median absolute deviations from the median library size
    rna.mult = (opts$mad_cutoff * mad(cellmd$rna.size))
    prot.mult = (opts$mad_cutoff * mad(cellmd$prot.size))
    rna.lower = median(cellmd$rna.size) - rna.mult
    rna.upper = median(cellmd$rna.size) + rna.mult
    prot.lower = median(cellmd$prot.size) - prot.mult
    prot.upper = median(cellmd$prot.size) + prot.mult

    # filter rows based on droplet qualty control metrics
    qc_cells = rownames(
        cellmd[cellmd$prot.size > prot.lower & 
                cellmd$prot.size < prot.upper & 
                cellmd$rna.size > rna.lower & 
                cellmd$rna.size < rna.upper & 
                cellmd$mt.prop < opts$max_mt_prop, ]
    )

    cell.adt.raw = as.matrix(prot[ , qc_cells])
    cell.rna.raw = rna[ ,qc_cells]
    cellmd = cellmd[qc_cells, ]

    # define isotype controls 
    isotype.controls = rownames(prot)[grep(opts$isotype_pattern, rownames(prot))]

    # normalize and denoise with dsb with 
    cells.dsb.norm = DSBNormalizeProtein(
        cell_protein_matrix = cell.adt.raw, 
        empty_drop_matrix = background.adt.mtx, 
        denoise.counts = TRUE, 
        use.isotype.control = TRUE, 
        isotype.control.name.vec = isotype.controls,
        return.stats = TRUE
    )

    # integrating with Seurat
    stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
    stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))

    # create Seurat object note: min.cells is a gene filter, not a cell filter
    s = Seurat::CreateSeuratObject(
            counts = cell.rna.raw, 
            meta.data = cellmd,
            assay = "RNA", 
            min.cells = 20
        )

    # add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot
    s[["CITE"]] = CreateAssayObject(data = cells.dsb.norm$dsb_normalized_matrix)

    # define proteins to use in clustering (non-isotype controls)
    prots = rownames(s@assays$CITE)[!grepl(opts$isotype_pattern, rownames(prot))]

    DefaultAssay(s) = 'CITE'

    # cluster and run umap 
    s = Seurat::FindNeighbors(
        object = s,
        dims = NULL,
        assay = 'CITE',
        features = prots,
        k.param = 30,
        verbose = TRUE
    )

    # direct graph clustering 
    s = Seurat::FindClusters(
        object = s,
        resolution = 1, 
        algorithm = 3, 
        graph.name = 'CITE_snn', 
        verbose = FALSE
    )

    # umap (optional)
    s = Seurat::RunUMAP(
        object = s,
        assay = "CITE",
        features = prots,
        seed.use = 1990,
        min.dist = 0.2,
        n.neighbors = 30,
        verbose = FALSE
    )

    # make results dataframe 
    d = cbind(
        s@meta.data, 
        as.data.frame(t(s@assays$CITE@data))
                # s@reductions$umap@cell.embeddings)
    )

    # calculate the median protein expression separately for each cluster 
    adt_plot = d %>% 
        dplyr::group_by(CITE_snn_res.1) %>% 
        dplyr::summarize_at(.vars = prots, .funs = median) %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames("CITE_snn_res.1")

    png(file = heatmap.png, height = 15, width = 8, res = 300, units = 'in')
    # plot a heatmap of the average dsb normalized values for each cluster
    pheatmap::pheatmap(t(adt_plot), 
                        color = viridis::viridis(25, option = "B"), 
                        fontsize_row = 8, border_color = NA)
    dev.off()

    log_message(sprintf("dsb: %d cells kept, %d background drops, %d isotype controls",
                        length(qc_cells), length(background_drops), length(isotype.controls)), "INFO")
    list(seurat = s, cells.dsb.norm = cells.dsb.norm, cellmd = cellmd,
         n_background = length(background_drops))
}

# -- MAIN ----------------------------------------------------------------------
res <- run_dsb(
  raw.path    = smk@params$raw_matrix,
  filt.path   = smk@params$filtered_matrix,
  opts        = smk@params$opts,
  mt.pattern  = smk@params$mt_pattern,
  scat.png    = smk@output$adt_scat_png,
  heatmap.png = smk@output$adt_heatmap_png
)
qs2::qs_save(object = res[c("seurat", "cells.dsb.norm", "cellmd")], file = smk@output$adt_dsb)
write.table(
  data.frame(sample = smk@wildcards$sample, dsb_cells = nrow(res$cellmd), background_drops = res$n_background),
  smk@output$summary, sep = "\t", quote = FALSE, row.names = FALSE
)
log_message("dsb done", "DONE")
