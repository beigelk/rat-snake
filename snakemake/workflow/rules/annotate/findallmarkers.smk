def findallmarkers_prefix(wc):
    return f"{PROJECT}__{wc.pop_name}__{wc.assay}_res.{wc.resolution}"

checkpoint findallmarkers:
    input:
        seurat_obj = compiled_seurat_input
    params:
        pop_name       = lambda wc: wc.pop_name,
        graph_name     = lambda wc: POPULATIONS[wc.pop_name]["graph_name"],
        assay          = lambda wc: wc.assay,
        resolution     = lambda wc: float(wc.resolution),
        resolution_str = lambda wc: wc.resolution, 
        prefix         = findallmarkers_prefix,
        r_library      = config['resources']['r_library'],
        smk_obj        = DEBUG_IDX['annotate/findallmarkers']['debug'] +
            "findallmarkers__{pop_name}__{assay}_res.{resolution}.rds"
    output:
        raw_tsv = pop_outdir("{pop_name}", "findallmarkers", "tables/raw") +
            f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__markers_raw.tsv"
    envmodules:
        "gcc",
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/findallmarkers']['log'] +
            "findallmarkers__{pop_name}__{assay}_res.{resolution}.log"
    script:
        RULE_IDX['annotate/findallmarkers']['r']
