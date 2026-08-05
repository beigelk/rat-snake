import json
rule metadata:
    input:
        seurat_obj      = harmonized_seurat
    params:
        metadata_json = json.dumps({
            d: config["datasets"][d]["metadata"]["input"]
            for d in config["datasets"]
        }),
        r_library = config['resources']['r_library'],
        smk_obj   = DEBUG_IDX['annotate/metadata']['debug'] + "metadata.rds"
    output:
        metadata_obj = metadata_obj
    envmodules:
        "R/4.5.3"
    log:
        logfile   = LOG_IDX['annotate/metadata']['log'] + "metadata.log"
    script:
        RULE_IDX['annotate/metadata']['r']

