COMPILE_TARGETS = list(POPULATIONS.keys())

rule compile:
    input:
        unpack(compile_inputs),
        metadata_obj = metadata_obj
    params:
        pop_name          = lambda wc: wc.pop_name,
        cluster_meta_json = lambda wc: json.dumps(population_cluster_metadata(wc.pop_name)),
        r_library         = config['resources']['r_library'],
        smk_obj           = DEBUG_IDX['annotate/compile']['debug'] +
            "compile__{pop_name}.rds"
    output:
        compiled_obj = pop_outdir("{pop_name}", "compile", "objects") +
            "{pop_name}_seurat.qs2"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/compile']['log'] + "compile__{pop_name}.log"
    script:
        RULE_IDX['annotate/compile']['r']

