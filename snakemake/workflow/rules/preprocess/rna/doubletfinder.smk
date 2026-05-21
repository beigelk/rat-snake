rule doubletfinder:
    input:
        soupx_object = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/objects/soupx__{sample}.qs2"
    params:
        r_library         = config['resources']['r_library'],
        sample_name       = "{sample}",
        min_feat_cutoff   = config['preprocess']['rna']['doubletfinder']['min_feat_cutoff'],
        max_feat_cutoff   = config['preprocess']['rna']['doubletfinder']['max_feat_cutoff'],
        percent_mt_cutoff = config['preprocess']['rna']['doubletfinder']['percent_mt_cutoff']
    output:
        doublet_id_file    = RES_IDX['preprocess/rna/doubletfinder']['outdir'] + "{dataset}/tables/{sample}_doublet_ids.txt",
        doublet_pkplot_pdf = RES_IDX['preprocess/rna/doubletfinder']['outdir'] + "{dataset}/figures/{sample}_doublet_pK.pdf",
        doublet_qc_png     = RES_IDX['preprocess/rna/doubletfinder']['outdir'] + "{dataset}/figures/{sample}_doublet_QC.png"
    log: 
        logfile = LOG_IDX['preprocess/rna/doubletfinder']['log'] + "{dataset}/{sample}.log"
    envmodules:
        "R/4.5.3"
    resources:
        tmpdir = "tmp/"
    script:
        RULE_IDX['preprocess/rna/doubletfinder']['r']