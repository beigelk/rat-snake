rule dsb:
    input:
        filt_path = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix",
        raw_path  = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/multi/count/raw_feature_bc_matrix"
    params:
        r_library   = config['resources']['r_library'],
        sample_name = "{sample}"
    output:
        adt_dsb         = RES_IDX['preprocess/adt/dsb']['outdir'] + "{dataset}/objects/{sample}_dsb.qs2",
        adt_scat_png    = RES_IDX['preprocess/adt/dsb']['outdir'] + "{dataset}/figures/{sample}_ADT_ScatterPlot.png",
        adt_heatmap_png = RES_IDX['preprocess/adt/dsb']['outdir'] + "{dataset}/figures/{sample}_ADT_Heatmap.png"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['preprocess/adt/dsb']['log'] + "{dataset}/{sample}.log"
    script:
        RULE_IDX['preprocess/adt/dsb']['r']