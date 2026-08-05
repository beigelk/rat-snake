wildcard_constraints:
    pop_name   = "|".join(POPULATIONS.keys()),
    graph_name = "|".join({p["graph_name"] for p in POPULATIONS.values()}),
    resolution = r"\d+\.?\d*",

rule cluster:
    input:
        seurat_obj = lambda wc: population_seurat_obj(wc.pop_name)
    params:
        run_config = lambda wc: POPULATIONS[wc.pop_name],
        pop_name   = lambda wc: wc.pop_name,
        graph_name = lambda wc: POPULATIONS[wc.pop_name]["graph_name"],
        resolution = lambda wc: float(wc.resolution),
        r_library  = config['resources']['r_library'],
        smk_obj    = DEBUG_IDX['annotate/cluster']['debug'] +
            "cluster__{pop_name}__{graph_name}_res.{resolution}.rds"
    output:
        results_tsv = pop_outdir("{pop_name}", "cluster", "tables") +
            "cluster__{pop_name}__{graph_name}_res.{resolution}.tsv"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/cluster']['log'] +
            "cluster__{pop_name}__{graph_name}_res.{resolution}.log"
    script:
        RULE_IDX['annotate/cluster']['r']