# ==============================================================================
# MODULE: annotate
# Cell metadata (sample sheets + clonotypes), then population by population:
# clusters at every resolution, subsets cut from a parent's clusters, one
# compiled object per population, markers and clustree.
# Populations come from `annotate.populations`; every path is built by
# workflow/lib/rat_snake/populations.py.
# ==============================================================================

import re

ANN_SCRIPTS = os.path.join(MODULES_DIR, "annotate", "scripts")
ANN_CFG = config["annotate"]
ANN_ENV = r_env("annotate")

def _alternatives(names):
    """Regex matching exactly one of the names (longest first)."""
    return "|".join(re.escape(n) for n in sorted(set(names), key=len, reverse=True)) or "NONE"


wildcard_constraints:
    population=_alternatives(POPS.names),
    res=r"\d+(\.\d+)?",
    assay=_alternatives(a for p in POPS.names for a in POPS.assays(p)),
    tset=_alternatives(t for p in POPS.names for a in POPS.assays(p) for t in POPS.threshold_sets(p, a)),


def dataset_metadata():
    """Each dataset's metadata block (name defaults to the dataset), or None."""
    out = {}
    for d in DATASETS:
        m = config["datasets"][d].get("metadata")
        out[d] = {**m, "name": m.get("name") or d} if m else None
    return out


rule metadata:
    """Build the cell metadata table: sample spreadsheets matched to samples, plus scRepertoire clonotypes."""
    input:
        seurat_obj=CTX.integrated_object,
        clones=[P.result("preprocess", "tcr_clones", d, "objects", f"{s}_combineTCR.qs2") for d, s in USABLE]
        if TCR_ENABLED
        else [],
    output:
        metadata_obj=CTX.metadata_table,
        tsv=P.result("annotate", "metadata", "cell_metadata.tsv"),
        matches=P.result("annotate", "metadata", "sample_matches.tsv"),
    params:
        metadata=dataset_metadata(),
        match_max_dist=ANN_CFG["metadata"]["match_max_dist"],
        tcr=config["preprocess"]["tcr"],
    log:
        P.log("annotate", "metadata"),
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "metadata.R")


rule cluster:
    """Cluster one population at one resolution (FindClusters on its WNN graph) and save the clusters."""
    input:
        seurat_obj=lambda wc: POPS.input_object(wc.population),
    output:
        results_tsv=POPS.cluster_tsv("{population}", "{res}"),
    params:
        pop_name=lambda wc: wc.population,
        graph_name=lambda wc: POPS.graph(wc.population),
        resolution=lambda wc: wc.res,
    log:
        P.log("annotate", "cluster", "{population}", "{res}"),
    benchmark:
        P.benchmark("annotate", "cluster", "{population}", "{res}")
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "cluster.R")


rule subset:
    """Cut a population from its parent's clusters at one resolution and re-run WNN and UMAP."""
    input:
        seurat_obj=lambda wc: POPS.input_object(POPS.subset_spec(wc.population)["base_name"]),
        cluster_results=lambda wc: POPS.cluster_tsv(
            POPS.subset_spec(wc.population)["base_name"], POPS.subset_spec(wc.population)["base_resolution"]
        ),
    output:
        seurat_sub_obj=POPS.subset_object("{population}"),
    params:
        pop_name=lambda wc: wc.population,
        graph_name=lambda wc: POPS.graph(wc.population),
        subset=lambda wc: POPS.subset_spec(wc.population),
        wnn=lambda wc: POPS.wnn(wc.population, ANN_CFG["subset_wnn"]),
    log:
        P.log("annotate", "subset", "{population}"),
    benchmark:
        P.benchmark("annotate", "subset", "{population}")
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "subset.R")


rule compile:
    """One object per population with every resolution's clusters and the cell metadata; cells per cluster."""
    input:
        target_obj=lambda wc: POPS.input_object(wc.population),
        cluster_results=lambda wc: [POPS.cluster_tsv(wc.population, r) for r in POPS.resolutions(wc.population)],
        metadata_obj=CTX.metadata_table,
    output:
        compiled_obj=POPS.compiled_object("{population}"),
        cluster_sizes=POPS.cluster_sizes("{population}"),
    params:
        pop_name=lambda wc: wc.population,
    log:
        P.log("annotate", "compile", "{population}"),
    benchmark:
        P.benchmark("annotate", "compile", "{population}")
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "compile.R")


rule findallmarkers:
    """Positive markers of every cluster (FindAllMarkers) for one population, resolution and assay."""
    input:
        seurat_obj=lambda wc: POPS.compiled_object(wc.population),
    output:
        raw_tsv=POPS.markers_raw("{population}", "{res}", "{assay}"),
    params:
        cluster_col=lambda wc: POPS.cluster_column(wc.population, wc.res),
        assay=lambda wc: wc.assay,
    retries: 2
    log:
        P.log("annotate", "findallmarkers", "{population}", "{res}", "{assay}"),
    benchmark:
        P.benchmark("annotate", "findallmarkers", "{population}", "{res}", "{assay}")
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "findallmarkers.R")


rule filtermarkers:
    """Keep the markers passing one threshold set (log2FC, adjusted p value, pct.1)."""
    input:
        raw_tsv=rules.findallmarkers.output.raw_tsv,
    output:
        filtered_tsv=POPS.markers_filtered("{population}", "{res}", "{assay}", "{tset}"),
    params:
        thresholds=lambda wc: POPS.threshold_sets(wc.population, wc.assay)[wc.tset],
    log:
        P.log("annotate", "filtermarkers", "{population}", "{res}", "{assay}", "{tset}"),
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "filtermarkers.R")


rule clustree:
    """clustree of a population's resolutions and a UMAP of its clusters at each one."""
    input:
        seurat_obj=lambda wc: POPS.compiled_object(wc.population),
    output:
        png=POPS.clustree("{population}"),
        umaps=directory(P.result("annotate", "populations", "{population}", "clustree", "umaps")),
    params:
        pop_name=lambda wc: wc.population,
        graph_name=lambda wc: POPS.graph(wc.population),
        project=PROJECT,
        reduction=config["visualization"]["reduction"],
    log:
        P.log("annotate", "clustree", "{population}"),
    conda:
        ANN_ENV.conda
    container:
        ANN_ENV.container
    script:
        os.path.join(ANN_SCRIPTS, "clustree.R")
