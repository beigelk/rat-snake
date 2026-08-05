rule subset:
    input:
        seurat_obj      = lambda wc: population_seurat_obj(POPULATIONS[wc.pop_name]["subset"]["base_name"]),
        cluster_results = lambda wc: cluster_tsv_path(
            POPULATIONS[wc.pop_name]["subset"]["base_name"],
            POPULATIONS[wc.pop_name]["subset"]["base_resolution"]
        )
    params:
        run_config      = lambda wc: POPULATIONS[wc.pop_name],
        pop_name        = lambda wc: wc.pop_name,
        graph_name      = lambda wc: POPULATIONS[wc.pop_name]["graph_name"],
        base_name       = lambda wc: POPULATIONS[wc.pop_name]["subset"]["base_name"],
        base_graph_name = lambda wc: POPULATIONS[POPULATIONS[wc.pop_name]["subset"]["base_name"]]["graph_name"],
        base_resolution = lambda wc: float(POPULATIONS[wc.pop_name]["subset"]["base_resolution"]),
        base_clusters   = lambda wc: parse_clusters_list(POPULATIONS[wc.pop_name]["subset"]["base_clusters"]),
        r_library       = config['resources']['r_library'],
        smk_obj         = DEBUG_IDX['annotate/subset']['debug'] +
            "subset__{pop_name}__{graph_name}.rds"
    output:
        seurat_sub_obj  = pop_outdir("{pop_name}", "subset", "objects") +
            "subset__{pop_name}__{graph_name}.qs2"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/subset']['log'] +
            "subset__{pop_name}__{graph_name}.log"
    script:
        RULE_IDX['annotate/subset']['r']
        