
# -- SNAKEMAKE VARIABLES ------------------- #
# input
filt_path = snakemake@input[['filt_path']]
raw_path  = snakemake@input[['raw_path']]

# params
r_library   = snakemake@params[['r_library']]
sample_name = snakemake@params[['sample_name']]

# output
adt_dsb         = snakemake@output[['adt_dsb']]
adt_scat_png    = snakemake@output[['adt_scat_png']]
adt_heatmap_png = snakemake@output[['adt_heatmap_png']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)

# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# -- LOAD PACKAGES ------------------------- #
# analysis
library(dsb)
library(Seurat)

# utility
library(qs2)

# plotting
library(tidyverse)
library(ggplot2)

# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 1500 * 1024^2)

# -- FUNCTIONS ------------------------------ #
# Normalizing and denoising CITE-seq (ADT) data
run_dsb = function(raw.path, filt.path)
{  
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
    mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below

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

    png(file = adt_scat_png, height = 6, width = 10, res = 300, units = 'in')
    print(
        ggplot(md, aes(x = log10(n.gene), y = prot.size)) +
        theme_bw() + 
        geom_bin2d(bins = 300) + 
        scale_fill_viridis_c(option = "C") + 
        facet_wrap(~drop.class)
    )
    dev.off()

    background_drops = rownames(
        md[ md$prot.size > 1.5 & 
            md$prot.size < 3 & 
            md$rna.size < 2.5, ]
    ) 
    background.adt.mtx = as.matrix(prot[ , background_drops])

    # calculate statistical thresholds for droplet filtering.
    cellmd = md[md$drop.class == 'cell', ]

    # filter drops with + / - 3 median absolute deviations from the median library size
    rna.mult = (3*mad(cellmd$rna.size))
    prot.mult = (3*mad(cellmd$prot.size))
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
                cellmd$mt.prop < 0.14, ]
    )

    length(qc_cells)

    cell.adt.raw = as.matrix(prot[ , qc_cells])
    cell.rna.raw = rna[ ,qc_cells]
    cellmd = cellmd[qc_cells, ]

    # define isotype controls 
    isotype.controls = rownames(prot)[grep('Isotype', rownames(prot))]

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
    prots = rownames(s@assays$CITE)[!grepl('Isotype', rownames(prot))]

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

    png(file = adt_heatmap_png, height = 15, width = 8, res = 300, units = 'in')
    # plot a heatmap of the average dsb normalized values for each cluster
    pheatmap::pheatmap(t(adt_plot), 
                        color = viridis::viridis(25, option = "B"), 
                        fontsize_row = 8, border_color = NA)
    dev.off()

    dsb_output = list(s, cells.dsb.norm, cellmd) %>% setNames(c('seurat', 'cells.dsb.norm', 'cellmd'))

    qs2::qs_save(
        object = dsb_output,
        file = adt_dsb
    )
  
}

# -- MAIN ----------------------------------------- #
run_dsb(raw.path = raw_path, filt.path = filt_path)
