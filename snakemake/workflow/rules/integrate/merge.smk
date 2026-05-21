rule merge:
    input:
        seurat_objects = complete_merge_inputs
    output:
        merged_seurat = RES_IDX['integrate/merge']['outdir'] + "objects/merged_seurat.qs2"
    params:
        r_library = config['resources']['r_library']
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['integrate/merge']['log'] + "merge_seurat.log"
    script:
        RULE_IDX['integrate/merge']['r']