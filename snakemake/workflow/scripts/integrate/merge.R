
# -- SNAKEMAKE VARIABLES ------------------- #
# input
seurat_objects = snakemake@input[['seurat_objects']]

# params
r_library = snakemake@params[['r_library']]

# output
merged_seurat = snakemake@output[['merged_seurat']]

# log
logfile = snakemake@log[1][[1]]

# -- SET LIBRARY PATH ----------------------- #
.libPaths(r_library)


# -- SET UP LOG ----------------------------- #
log <- file(logfile, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# -- LOG INPUTS ------------------------------ #
print(snakemake)
dput(seurat_objects)

# -- LOAD PACKAGES --------------------------- #
# analysis
library(Seurat)

# utility
library(qs2)

# plot
library(tidyverse)
library(ggplot2)

# -- ENVIRONMENT SETTINGS ------------------- #
set.seed(42)
options(future.globals.maxSize = 200000 * 1024^2)


# -- FUNCTIONS ------------------------------ #
make_seurat_list = function(seurat.file.list) {
    # seurat.obj.list: list of Seurat objects to merge
    seurat.obj.list = vector(mode = "list")

    for (filename in seurat.file.list){
        
        # Load the Seurat object
        seu = qs2::qs_read(file = filename)
        print(seu@project.name)
        print(seu)
        
        # Add to list
        seurat.obj.list[[seu@project.name]] = seu
        
        rm(seu)
    }

    # Removing low-quality samples
    ###################################################################
    # -----------------------------------------------------------------
    # samples.to.rm = read.table(
    #     file = exclusion_metadata, header = TRUE, sep = ','
    # ) %>%
    # filter(Reason == 'low_qual') %>%
    # pull(Sample)

    # # Remove spotty samples
    # seurat.obj.list = seurat.obj.list[!names(seurat.obj.list) %in% samples.to.rm]
    # -----------------------------------------------------------------
    ###################################################################

    return(seurat.obj.list)
}


make_merged_seurat = function(seurat.obj.list) {  
    # Merge the list into a merge object (Seurat 5)
    merged= merge(
        x = seurat.obj.list[[1]],
        y = seurat.obj.list[2:length(seurat.obj.list)],
        merge.data = TRUE,
        add.cell.ids = names(seurat.obj.list)
    )
  
  return(merged)

}

###########################################################################

seurat_obj_list = make_seurat_list(seurat_objects)
merged = make_merged_seurat(seurat_obj_list)

qs2::qs_save(
    object = merged,
    file = merged_seurat
)
