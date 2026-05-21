
vdj_contig_anno = ''
combined_tcr = ''
outfile_base = ''

#########
# Analysis
library(Seurat)
library(scRepertoire)

# Utility
library(qs2)
library(tidyverse)
library(ggplot2)
#########

#########
set.seed(42)
#########

source('Rscripts/utils.R', print.eval = FALSE)

# Create the figure dir from the outfile base
dir.create(dirname(outfile_base), showWarnings = FALSE, recursive = TRUE)

# ###############################################################################################################################
# FUNCTIONS
############################################################################

# LOAD ALL TCR DATA (cellranger outs) INTO A LIST
# --------------------------------------------------------------------------
# This works for files structured as:
# seq.dir
#    ├── sample_1
#    │      └── outs/per_sample_outs/sample/vdj_t/filtered_contig_annotations.csv
#    ├── sample_2
#    │      └── outs/per_sample_outs/sample/vdj_t/filtered_contig_annotations.csv
#    └── sample_3
#           └── outs/per_sample_outs/sample/vdj_t/filtered_contig_annotations.csv

load_tcr = function(seq.dir) {

    # From the main directory, get a list of all the sample folders
    sample.list = list.files(seq.dir)

    # Init a list
    contig.anno = vector(mode = "list")

    # For each sample,
    # load the TCR data file 'filtered_contig_annotations.csv'
    # and add to the contig.anno list
    for (sample in sample.list) {
    
    # Get the path to TCR data file for the sample
    vdj.file = paste0(seq.dir, sample, '/outs/per_sample_outs/', sample, '/vdj_t/filtered_contig_annotations.csv')
    
    # Read in CSV of the TCR data for the sample
    sample.tcr = read.csv(vdj.file) %>%
        mutate(sample = sample)
    
    # Add to list of TCR data
    contig.anno[[sample]] = sample.tcr

    }

    return(contig.anno)

}
# --------------------------------------------------------------------------

# UPDATE BARCODE NAMES TO MATCH THE TCR DATA ('BARCODE_SAMPLENAME')
# --------------------------------------------------------------------------
change_barcodes = function(seu.obj, sample.ids) {

    # seu.obj: Seurat object
    # sample.ids: the metadata column in seu.obj that has sample names 

    # Make a new column in the Seurat obj metadata of the original barcodes
    seu.obj = AddMetaData(seu.obj, metadata = colnames(seu.obj), col.name = 'orig_barcode')

    # scRepertoire adds the sample name as a PREFIX to the barcode,
    # but Seurat adds a number as a suffux.
    # Here, we will gsub out the number suffix from the Seurat barcode
    # and add the sample name as a prefix.
    new.bcs = seu.obj@meta.data %>%
              select(!!sym(sample.ids), orig_barcode) %>%
              mutate(new_barcode = paste0(orig.ident, '_', gsub('_.*$', '', orig_barcode)))

    colnames(seu.obj) = new.bcs$new_barcode

    return(seu.obj)

}
# --------------------------------------------------------------------------
############################################################################


# MAIN
# --------------------------------------------------------------------------

# Load the TCR data
# tcrseq_dir should be path to 10x data files output from CellRanger multi
contig_anno = load_tcr(tcrseq_dir)

# Combine TCR data in scRepertoire
combined_TCR = combineTCR(
    contig_anno,
    samples = names(contig_anno),
    removeNA = FALSE, 
    removeMulti = FALSE, 
    filterMulti = FALSE
)

# Change barcodes to match scRepertoire convention
seurat = change_barcodes(seurat, 'orig.ident')

# Add the TCR data to the Seurat object
seurat = combineExpression(combined_TCR, seurat)

DimPlot(seurat, group.by = 'cloneSize')

qs2::qs_save(
    combined_TCR,
    file = screpertoire_outfile
)

qs2::qs_save(
    seurat,
    file = seurat_outfile
)
# --------------------------------------------------------------------------
