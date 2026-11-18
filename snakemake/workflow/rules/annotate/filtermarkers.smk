def filtermarkers_prefix(wc):
    return f"{PROJECT}__{wc.pop_name}__{wc.assay}_res.{wc.resolution}_{wc.threshold_set}"

checkpoint filtermarkers:
    input:
        raw_tsv = lambda wc: markers_raw_tsv_path(wc.pop_name, wc.resolution, wc.assay)
    params:
        pop_name       = lambda wc: wc.pop_name,
        assay          = lambda wc: wc.assay,
        threshold_set  = lambda wc: wc.threshold_set,
        lfc_threshold  = lambda wc: pop_threshold_values(wc.pop_name, wc.assay, wc.threshold_set)['lfc_threshold'],
        padj_threshold = lambda wc: pop_threshold_values(wc.pop_name, wc.assay, wc.threshold_set)['padj_threshold'],
        pct1_threshold = lambda wc: pop_threshold_values(wc.pop_name, wc.assay, wc.threshold_set)['pct1_threshold'],
        prefix         = filtermarkers_prefix,
        r_library      = config['resources']['r_library'],
        smk_obj        = DEBUG_IDX['annotate/findallmarkers']['debug'] +
            "filtermarkers__{pop_name}__{assay}_res.{resolution}__{threshold_set}.rds"
    output:
        filtered_tsv = pop_outdir("{pop_name}", "findallmarkers", "tables/filtered") +
            f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{threshold_set}}__markers_filtered.tsv"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/findallmarkers']['log'] +
            "filtermarkers__{pop_name}__{assay}_res.{resolution}__{threshold_set}.log"
    script:
        RULE_IDX['annotate/filtermarkers']['r']