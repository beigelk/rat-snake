rule construct:
    input:
        soupx_object    = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/objects/soupx__{sample}.qs2",
        doublet_id_file = RES_IDX['preprocess/rna/doubletfinder']['outdir'] + "{dataset}/tables/{sample}_doublet_ids.txt",
        adt_dsb_object  = RES_IDX['preprocess/adt/dsb']['outdir'] + "{dataset}/objects/{sample}_dsb.qs2",
        # vdj_contig_anno = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv"
    params:
        r_library             = config['resources']['r_library'],
        sample_name           = '{sample}',
        rna_min_feat_cutoff   = config['preprocess']['multimodal']['construct']['rna']['rna_min_feat_cutoff'],
        rna_max_feat_cutoff   = config['preprocess']['multimodal']['construct']['rna']['rna_max_feat_cutoff'],
        rna_percent_mt_cutoff = config['preprocess']['multimodal']['construct']['rna']['rna_percent_mt_cutoff'],
        adt_feat_min          = config['preprocess']['multimodal']['construct']['adt']['adt_feat_min'],
        adt_isocount_max      = config['preprocess']['multimodal']['construct']['adt']['adt_isocount_max'],
    output:
        seurat_object        = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/objects/{sample}_Seurat.qs2",
        rna_filt_cutoffs     = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/tables/{sample}_RNA_filtering_cutoffs.txt",
        rna_sig_pcs          = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/tables/{sample}_RNA_sig_pcs.txt",
        rna_qcplots_prefilt  = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Prefilt_QCPlots.pdf",
        rna_umap_prefilt     = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Prefilt_UMAPPlot.png",
        rna_vln_prefilt      = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Prefilt_VlnPlot.png",
        rna_scat_prefilt     = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Prefilt_FeatScatter.png",
        rna_qcplots_postfilt = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Postfilt_QCPlots.pdf",
        rna_umap_postfilt    = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Postfilt_UMAPPlot.png",
        rna_vln_postfilt     = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Postfilt_VlnPlot.png",
        rna_scat_postfilt    = RES_IDX['preprocess/multimodal/construct']['outdir'] + "{dataset}/figures/{sample}_Postfilt_FeatScatter.png"
    envmodules:
        "R/4.5.3"
    log: 
        logfile = LOG_IDX['preprocess/multimodal/construct']['log'] + "{dataset}/{sample}.log"
    resources:
        tmpdir = 'tmp/'
    script:
        RULE_IDX['preprocess/multimodal/construct']['r']