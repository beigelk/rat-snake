rule tcr:
    input:
        vdj_contig_anno = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv"
    params:
        sample_name = "{sample}"
    output:
        adt_dsb         = RES_IDX['preprocess/tcr/clones']['outdir'] + "objects/{sample}_dsb.qs2",
        adt_scat_png    = RES_IDX['preprocess/tcr/clones']['outdir'] + "figures/{sample}_ADT_ScatterPlot.png",
        adt_heatmap_png = RES_IDX['preprocess/tcr/clones']['outdir'] + "figures/{sample}_ADT_Heatmap.png"
    log:
        logfile = LOG_IDX['preprocess/tcr/clones']['log'] + "{sample}.log"
    script:
        RULE_IDX['preprocess/tcr/clones']['r']