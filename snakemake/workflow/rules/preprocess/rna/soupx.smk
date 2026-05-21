rule soupx:
    input:
        filt_path = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix",
        raw_path  = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/multi/count/raw_feature_bc_matrix"
    params:
        r_library   = config['resources']['r_library'],
        sample_name = "{sample}",
        npcs        = config['preprocess']['rna']['soupx']['npcs'],
        resolution  = config['preprocess']['rna']['soupx']['resolution']
    output:
        soupx_object  = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/objects/soupx__{sample}.qs2",
        soupx_figure  = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/figures/soupx__{sample}.pdf"
    envmodules:
        "R/4.5.3"
    log: 
        logfile = LOG_IDX['preprocess/rna/soupx']['log'] + "{dataset}/{sample}.log"
    script:
        RULE_IDX['preprocess/rna/soupx']['r']
