# ==============================================================================
# MODULE: report
# One report of the run as HTML, Word and a plain Word file that pastes
# cleanly into Notion (snake_core.report): samples, Cell Ranger metrics,
# preprocessing QC, TCR clonotypes, the integrated object's metadata plots and,
# per population, cluster sizes, clustree, top markers, dot plots and marker
# plots. Full tables stay in results/.
# ==============================================================================

REPORT_SCRIPTS = os.path.join(MODULES_DIR, "report", "scripts")
REPORT_ENV = module_env("report", config, MODULES_DIR)
REPORT_CFG = config["report"]
REPORT_FILES = {"html": "report.html", "docx": "report.docx", "notion": "report_notion.docx"}


def report_resolutions(pop):
    """Resolutions of a population shown in the report (report.resolutions, default all)."""
    chosen = (REPORT_CFG.get("resolutions") or {}).get(pop)
    allres = POPS.resolutions(pop)
    if not chosen:
        return allres
    chosen = [round(float(r), 4) for r in chosen]
    unknown = [r for r in chosen if r not in allres]
    if unknown:
        raise WorkflowError(f"report.resolutions.{pop}: {unknown} are not resolutions of '{pop}'")
    return chosen


def report_plan():
    """Per population: resolutions shown, and the marker tables and plots of each."""
    tset = REPORT_CFG["threshold_set"]
    plan = {}
    for pop in POPS.names:
        res_list = report_resolutions(pop)
        plan[pop] = {
            "resolutions": res_list,
            "cluster_sizes": POPS.cluster_sizes(pop),
            "clustree": POPS.clustree(pop),
            "markers": {
                str(r): {a: POPS.markers_filtered(pop, r, a, tset) for a in POPS.assays(pop) if tset in POPS.threshold_sets(pop, a)}
                for r in res_list
            },
            "plots": {
                str(r): [
                    {"assay": a, "marker_set": m, "dotplot": POPS.dotplot(pop, r, a, m), "markerplot": POPS.markerplot(pop, r, a, m)}
                    for a in POPS.assays(pop)
                    for m in POPS.marker_sets(pop, a)
                ]
                for r in res_list
            },
        }
    return plan


REPORT_PLAN = report_plan()


def _files(obj):
    """Every path in the plan (strings ending in a file extension)."""
    if isinstance(obj, dict):
        return [f for v in obj.values() for f in _files(v)]
    if isinstance(obj, list):
        return [f for v in obj for f in _files(v)]
    return [obj] if isinstance(obj, str) and os.sep in obj else []


def per_sample(rule_name, kind, suffix):
    """One file per usable sample under results/preprocess/{rule_name}/{dataset}/{kind}/."""
    return [P.result("preprocess", rule_name, d, kind, f"{s}{suffix}") for d, s in USABLE]


rule build_report:
    """Write the run's report as Markdown and the formats in report.formats (HTML, DOCX, Notion DOCX)."""
    input:
        status=P.result("setup", "sample_status.tsv"),
        cellranger_metrics=P.result("cellranger", "metrics", "cellranger_metrics.tsv"),
        preprocess_summary=P.result("preprocess", "preprocess_summary.tsv"),
        doublet_qc=per_sample("doubletfinder", "figures", "_doublet_QC.png"),
        adt_scatter=per_sample("dsb", "figures", "_ADT_ScatterPlot.png"),
        postfilt_vln=per_sample("construct", "figures", "_Postfilt_VlnPlot.png"),
        clone_sizes=per_sample("tcr_clones", "figures", "_clone_sizes.png") if TCR_ENABLED else [],
        sample_matches=P.result("annotate", "metadata", "sample_matches.tsv"),
        metadata_plots=P.result("visualization", "metadata", "metadata_plots.tsv"),
        populations=_files(REPORT_PLAN),
    output:
        md=P.report("report.md"),
        formats=[P.report(REPORT_FILES[f]) for f in REPORT_CFG["formats"]],
    params:
        project=config["project"],
        options=REPORT_CFG,
        usable=[list(x) for x in USABLE],
        modalities=MODALITIES,
        tcr=TCR_ENABLED,
        plan=REPORT_PLAN,
        settings={
            "preprocess": config["preprocess"],
            "integrate": config["integrate"],
            "visualization": config["visualization"],
        },
    log:
        P.log("report", "build_report"),
    conda:
        REPORT_ENV.conda
    container:
        REPORT_ENV.container
    script:
        os.path.join(REPORT_SCRIPTS, "build_report.py")
