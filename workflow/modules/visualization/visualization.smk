# ==============================================================================
# MODULE: visualization
# UMAPs and QC plots of the integrated object coloured by sample metadata, and
# per population, resolution, assay and marker set: dot plots and a combined
# UMAP / violin / feature plot. Marker sets come from each population's
# `plots.marker_sets` (topn = top markers of a filtered table, genelist = a
# file with one feature per line).
# ==============================================================================

VIS_SCRIPTS = os.path.join(MODULES_DIR, "visualization", "scripts")
VIS_CFG = config["visualization"]
VIS_ENV = r_env("visualization")


wildcard_constraints:
    mset=_alternatives(m for p in POPS.names for a in POPS.assays(p) for m in POPS.marker_sets(p, a)),


def marker_set_spec(wc):
    """The marker set's spec ({type: topn, threshold_set, top_n} or {type: genelist, file})."""
    return POPS.marker_sets(wc.population, wc.assay)[wc.mset]


rule metadata_plots:
    """UMAPs and QC violin/ridge plots of the integrated object (annotate.metadata.plots)."""
    input:
        seurat_obj=CTX.integrated_object,
        metadata_obj=CTX.metadata_table,
    output:
        figures=directory(P.result("visualization", "metadata", "figures")),
        index=P.result("visualization", "metadata", "metadata_plots.tsv"),
    params:
        plot_config=ANN_CFG["metadata"]["plots"],
        reduction=VIS_CFG["reduction"],
        project=PROJECT,
    log:
        P.log("visualization", "metadata_plots"),
    conda:
        VIS_ENV.conda
    container:
        VIS_ENV.container
    script:
        os.path.join(VIS_SCRIPTS, "metadata_plots.R")


rule dotplot:
    """Dot plot and clustered dot plot of one marker set at one resolution."""
    input:
        seurat_obj=lambda wc: POPS.compiled_object(wc.population),
        features_file=lambda wc: POPS.plot_features_input(wc.population, wc.res, wc.assay, wc.mset),
    output:
        png=POPS.dotplot("{population}", "{res}", "{assay}", "{mset}"),
    params:
        spec=marker_set_spec,
        cluster_col=lambda wc: POPS.cluster_column(wc.population, wc.res),
        prefix=lambda wc: POPS.plot_prefix(wc.population, wc.res, wc.assay, wc.mset),
        assay=lambda wc: wc.assay,
    log:
        P.log("visualization", "dotplot", "{population}", "{res}", "{assay}", "{mset}"),
    conda:
        VIS_ENV.conda
    container:
        VIS_ENV.container
    script:
        os.path.join(VIS_SCRIPTS, "dotplot.R")


rule markerplots:
    """Clusters UMAP, stacked violins and feature UMAPs of one marker set (one combined figure)."""
    input:
        seurat_obj=lambda wc: POPS.compiled_object(wc.population),
        features_file=lambda wc: POPS.plot_features_input(wc.population, wc.res, wc.assay, wc.mset),
    output:
        png=POPS.markerplot("{population}", "{res}", "{assay}", "{mset}"),
    params:
        spec=marker_set_spec,
        cluster_col=lambda wc: POPS.cluster_column(wc.population, wc.res),
        prefix=lambda wc: POPS.plot_prefix(wc.population, wc.res, wc.assay, wc.mset),
        assay=lambda wc: wc.assay,
        reduction=VIS_CFG["reduction"],
    log:
        P.log("visualization", "markerplots", "{population}", "{res}", "{assay}", "{mset}"),
    conda:
        VIS_ENV.conda
    container:
        VIS_ENV.container
    script:
        os.path.join(VIS_SCRIPTS, "markerplots.R")
