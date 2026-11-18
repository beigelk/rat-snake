def plot_input(wc):
    spec = pop_marker_set_spec(wc.pop_name, wc.assay, wc.marker_set)
    inputs = {"seurat_obj": compiled_seurat_input(wc)}
    if spec["type"] == "topn":
        inputs["features_file"] = markers_filtered_tsv_path(
            wc.pop_name, wc.resolution, wc.assay, spec["threshold_set"]
        )
    elif spec["type"] == "genelist":
        inputs["features_file"] = spec["file"]
    else:
        raise ValueError(f"Unknown marker_set type: {spec['type']}")
    return inputs

def plot_prefix(wc):
    return f"{PROJECT}__{wc.pop_name}__{wc.assay}_res.{wc.resolution}_{wc.marker_set}"

rule dotplot:
    input:
        unpack(plot_input)
    params:
        pop_name    = lambda wc: wc.pop_name,
        graph_name  = lambda wc: POPULATIONS[wc.pop_name]["graph_name"],
        assay       = lambda wc: wc.assay,
        resolution  = lambda wc: float(wc.resolution),
        marker_set  = lambda wc: wc.marker_set,
        spec        = lambda wc: pop_marker_set_spec(wc.pop_name, wc.assay, wc.marker_set),
        outfile_base = pop_outdir("{pop_name}", "plots", "dotplot"),
        prefix       = plot_prefix,
        r_library    = config['resources']['r_library'],
        smk_obj      = DEBUG_IDX['visualization/dotplot']['debug'] +
             f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__dotplot.rds"
    output:
        done = pop_outdir("{pop_name}", "plots", "dotplot") +
             f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__dotplot.done"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['visualization/dotplot']['log'] +
             f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__dotplot.log"
    script:
        RULE_IDX['visualization/dotplot']['r']

rule markerplots:
    input:
        unpack(plot_input)
    params:
        pop_name    = lambda wc: wc.pop_name,
        graph_name  = lambda wc: wc.graph_name,
        resolution  = lambda wc: float(wc.resolution),
        assay       = lambda wc: wc.assay,
        marker_set  = lambda wc: wc.marker_set,
        spec        = lambda wc: pop_marker_set_spec(wc.pop_name, wc.assay, wc.marker_set),
        outfile_base = pop_outdir("{pop_name}", "plots", "markerplots"),
        prefix       = plot_prefix,
        r_library    = config['resources']['r_library'],
        smk_obj      = DEBUG_IDX['visualization/markerplots']['debug'] +
            f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__umap+violin+feature.rds"
    output:
        done = pop_outdir("{pop_name}", "plots", "markerplots") +
            f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__umap+violin+feature.done"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['visualization/markerplots']['log'] +
            f"{PROJECT}__{{pop_name}}__{{assay}}_res.{{resolution}}__{{marker_set}}__umap+violin+feature.log"
    script:
        RULE_IDX['visualization/markerplots']['r']