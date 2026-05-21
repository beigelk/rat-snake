rule harmony:
    input:
        merged_seurat = merged_seurat
    output:
        harmonized_seurat = harmonized_seurat
    params:
        r_library = config['resources']['r_library'],
        adt_markers_file = config['integrate']['harmony']['adt_markers_file']
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['integrate/harmony']['log'] + "harmony_seurat.log"
    script:
        RULE_IDX['integrate/harmony']['r']