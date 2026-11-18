# ==============================================================================
# MODULE: preprocess
# Per sample: ambient RNA (SoupX), doublets (DoubletFinder), ADT
# normalisation (dsb), TCR clonotypes (scRepertoire) and the filtered RNA +
# ADT Seurat object (construct). Then one QC table across samples.
# ==============================================================================

PRE_SCRIPTS = os.path.join(MODULES_DIR, "preprocess", "scripts")
PRE_CFG = config["preprocess"]
PRE_ENV = r_env("preprocess")
DOUBLET_ENV = r_env("preprocess", "doublet")


def pre(rule_name, *parts):
    """Path under results/preprocess/{rule_name}/{dataset}/."""
    return P.result("preprocess", rule_name, "{dataset}", *parts)


rule soupx:
    """Remove ambient RNA with SoupX and build the sample's first Seurat object (RNA + ADT)."""
    input:
        cellranger=cellranger_input,
    output:
        soupx_object=pre("soupx", "objects", "soupx__{sample}.qs2"),
        soupx_figure=pre("soupx", "figures", "soupx__{sample}.pdf"),
        summary=pre("soupx", "tables", "{sample}_soupx_summary.tsv"),
    params:
        filtered_matrix=lambda wc: CTX.cellranger_file(wc.dataset, wc.sample, "filtered_matrix"),
        raw_matrix=lambda wc: CTX.cellranger_file(wc.dataset, wc.sample, "raw_matrix"),
        npcs=PRE_CFG["soupx"]["npcs"],
        resolution=PRE_CFG["soupx"]["resolution"],
    resources:
        soupx_jobs=1,
    log:
        P.log("preprocess", "soupx", "{dataset}", "{sample}"),
    benchmark:
        P.benchmark("preprocess", "soupx", "{dataset}", "{sample}")
    conda:
        PRE_ENV.conda
    container:
        PRE_ENV.container
    script:
        os.path.join(PRE_SCRIPTS, "soupx.R")


rule doubletfinder:
    """Pre-filter cells and call doublets with DoubletFinder (either of two runs calling a doublet counts)."""
    input:
        soupx_object=rules.soupx.output.soupx_object,
    output:
        doublet_id_file=pre("doubletfinder", "tables", "{sample}_doublet_ids.txt"),
        doublet_pkplot_pdf=pre("doubletfinder", "figures", "{sample}_doublet_pK.pdf"),
        doublet_qc_png=pre("doubletfinder", "figures", "{sample}_doublet_QC.png"),
        summary=pre("doubletfinder", "tables", "{sample}_doublet_summary.tsv"),
    params:
        **PRE_CFG["doubletfinder"],
        mt_pattern=PRE_CFG["mt_pattern"],
    resources:
        doubletfinder_jobs=1,
        tmpdir="tmp/",
    log:
        P.log("preprocess", "doubletfinder", "{dataset}", "{sample}"),
    benchmark:
        P.benchmark("preprocess", "doubletfinder", "{dataset}", "{sample}")
    conda:
        DOUBLET_ENV.conda
    container:
        DOUBLET_ENV.container
    script:
        os.path.join(PRE_SCRIPTS, "doubletfinder.R")


rule dsb:
    """Normalise and denoise ADT counts with dsb, using empty drops as background and isotype controls."""
    input:
        cellranger=cellranger_input,
    output:
        adt_dsb=pre("dsb", "objects", "{sample}_dsb.qs2"),
        adt_scat_png=pre("dsb", "figures", "{sample}_ADT_ScatterPlot.png"),
        adt_heatmap_png=pre("dsb", "figures", "{sample}_ADT_Heatmap.png"),
        summary=pre("dsb", "tables", "{sample}_dsb_summary.tsv"),
    params:
        filtered_matrix=lambda wc: CTX.cellranger_file(wc.dataset, wc.sample, "filtered_matrix"),
        raw_matrix=lambda wc: CTX.cellranger_file(wc.dataset, wc.sample, "raw_matrix"),
        opts=PRE_CFG["dsb"],
        mt_pattern=PRE_CFG["mt_pattern"],
    resources:
        dsb_jobs=1,
    log:
        P.log("preprocess", "dsb", "{dataset}", "{sample}"),
    benchmark:
        P.benchmark("preprocess", "dsb", "{dataset}", "{sample}")
    conda:
        PRE_ENV.conda
    container:
        PRE_ENV.container
    script:
        os.path.join(PRE_SCRIPTS, "dsb.R")


rule tcr_clones:
    """Build the sample's TCR clonotypes from Cell Ranger's filtered contigs with scRepertoire combineTCR."""
    input:
        cellranger=cellranger_input,
    output:
        clones=pre("tcr_clones", "objects", "{sample}_combineTCR.qs2"),
        table=pre("tcr_clones", "tables", "{sample}_clones.tsv"),
        figure=pre("tcr_clones", "figures", "{sample}_clone_sizes.png"),
        summary=pre("tcr_clones", "tables", "{sample}_tcr_summary.tsv"),
    params:
        contigs=lambda wc: CTX.cellranger_file(wc.dataset, wc.sample, "contigs"),
        opts=PRE_CFG["tcr"],
    log:
        P.log("preprocess", "tcr_clones", "{dataset}", "{sample}"),
    conda:
        PRE_ENV.conda
    container:
        PRE_ENV.container
    script:
        os.path.join(PRE_SCRIPTS, "tcr_clones.R")


rule construct:
    """Build the sample's filtered RNA + ADT Seurat object: RNA QC, doublet removal, dsb values, ADT QC."""
    input:
        soupx_object=rules.soupx.output.soupx_object,
        doublet_id_file=rules.doubletfinder.output.doublet_id_file,
        adt_dsb_object=rules.dsb.output.adt_dsb,
    output:
        seurat_object=pre("construct", "objects", "{sample}_Seurat.qs2"),
        rna_filt_cutoffs=pre("construct", "tables", "{sample}_RNA_filtering_cutoffs.txt"),
        rna_sig_pcs=pre("construct", "tables", "{sample}_RNA_sig_pcs.txt"),
        summary=pre("construct", "tables", "{sample}_construct_summary.tsv"),
        rna_qcplots_prefilt=pre("construct", "figures", "{sample}_Prefilt_QCPlots.pdf"),
        rna_umap_prefilt=pre("construct", "figures", "{sample}_Prefilt_UMAPPlot.png"),
        rna_vln_prefilt=pre("construct", "figures", "{sample}_Prefilt_VlnPlot.png"),
        rna_scat_prefilt=pre("construct", "figures", "{sample}_Prefilt_FeatScatter.png"),
        rna_qcplots_postfilt=pre("construct", "figures", "{sample}_Postfilt_QCPlots.pdf"),
        rna_umap_postfilt=pre("construct", "figures", "{sample}_Postfilt_UMAPPlot.png"),
        rna_vln_postfilt=pre("construct", "figures", "{sample}_Postfilt_VlnPlot.png"),
        rna_scat_postfilt=pre("construct", "figures", "{sample}_Postfilt_FeatScatter.png"),
    params:
        rna=PRE_CFG["construct"]["rna"],
        adt=PRE_CFG["construct"]["adt"],
        npcs=PRE_CFG["construct"]["npcs"],
        cluster_resolution=PRE_CFG["construct"]["cluster_resolution"],
        mt_pattern=PRE_CFG["mt_pattern"],
        isotype_pattern=PRE_CFG["dsb"]["isotype_pattern"],
    resources:
        construct_jobs=1,
        tmpdir="tmp/",
    log:
        P.log("preprocess", "construct", "{dataset}", "{sample}"),
    benchmark:
        P.benchmark("preprocess", "construct", "{dataset}", "{sample}")
    conda:
        PRE_ENV.conda
    container:
        PRE_ENV.container
    script:
        os.path.join(PRE_SCRIPTS, "construct.R")


def preprocess_tables(kind):
    """Every usable sample's summary table of one step."""
    files = {
        "soupx": ("soupx", "_soupx_summary.tsv"),
        "doubletfinder": ("doubletfinder", "_doublet_summary.tsv"),
        "dsb": ("dsb", "_dsb_summary.tsv"),
        "construct": ("construct", "_construct_summary.tsv"),
        "tcr": ("tcr_clones", "_tcr_summary.tsv"),
    }
    rule_name, suffix = files[kind]
    return [P.result("preprocess", rule_name, d, "tables", f"{s}{suffix}") for d, s in USABLE]


rule preprocess_summary:
    """Join every sample's SoupX, DoubletFinder, dsb, construct and TCR summaries into one QC table."""
    input:
        soupx=preprocess_tables("soupx"),
        doublets=preprocess_tables("doubletfinder"),
        dsb=preprocess_tables("dsb"),
        construct=preprocess_tables("construct"),
        tcr=preprocess_tables("tcr") if TCR_ENABLED else [],
    output:
        tsv=P.result("preprocess", "preprocess_summary.tsv"),
    log:
        P.log("preprocess", "preprocess_summary"),
    run:
        import pandas as pd

        def read_all(files):
            # Files are in USABLE order; the dataset column keeps equal sample names apart
            parts = [pd.read_csv(f, sep="\t").assign(dataset=d) for f, (d, _) in zip(files, USABLE)]
            return pd.concat(parts, ignore_index=True)

        table = read_all(input.construct)
        for part in ("soupx", "doublets", "dsb", "tcr"):
            files = input[part]
            if files:
                extra = read_all(files).drop(columns=["cells"], errors="ignore")
                table = table.merge(extra, on=["dataset", "sample"], how="left")
        cols = ["dataset", "sample"] + [c for c in table.columns if c not in ("dataset", "sample")]
        table[cols].to_csv(output.tsv, sep="\t", index=False)
