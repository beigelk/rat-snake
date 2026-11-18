"""build_report: the run's report of figures and tables.

Snakemake script (rule ``build_report``). Reads the tables and figures of the
setup, cellranger, preprocess, annotate and visualization modules and writes
report.md plus the formats in ``report.formats`` (report.html, report.docx,
report_notion.docx) with snake_core.report.

Sections:
  1 Summary           samples, cells, populations
  2 Samples           status of every sample in the sample tables
  3 Cell Ranger       key metrics per sample
  4 Preprocessing     cells kept at each step, ambient RNA, doublets, ADT
                      background; per-sample QC figures
  5 TCR clonotypes    clonotypes per sample and clone size figures
  6 Integrated object sample metadata matches and metadata UMAPs / QC plots
  7 Populations       per population: cluster sizes, clustree, and per shown
                      resolution the top markers, dot plots and marker plots
  8 Settings          the main analysis options
Params:  project, options (``report``), usable, modalities, tcr, plan
         (populations, see report.smk), settings
"""

import os

import pandas as pd
import yaml
from snake_core.logging import log_message, redirect_to_log
from snake_core.report import Report

redirect_to_log(snakemake.log[0])
inp = snakemake.input
par = snakemake.params
opts = dict(par.options or {})
max_rows = opts.get("max_table_rows", 40)
top_n = opts.get("top_markers", 5)
project = par.project


def read(path, **kw):
    """Read a TSV, or an empty table when it is missing or empty."""
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def sample_of(path, suffix):
    """``dataset/sample`` of a per-sample file .../{dataset}/{kind}/{sample}{suffix}."""
    dataset = os.path.basename(os.path.dirname(os.path.dirname(path)))
    return f"{dataset}/{os.path.basename(path)[: -len(suffix)]}"


def top_markers(path, n):
    """Top n markers per cluster of a filtered marker table (in its order)."""
    df = read(path)
    if df.empty:
        return df
    keep = [c for c in ["cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"] if c in df.columns]
    top = df.groupby("cluster", sort=False).head(n)[keep].copy()
    for col in ("p_val_adj",):
        if col in top.columns:
            top[col] = top[col].map(lambda v: f"{v:.2e}")
    return top


status = read(inp.status)
qc = read(inp.preprocess_summary)
plan = par.plan

r = Report(
    opts.get("title") or f"{project['name']} {project.get('version', '')}".strip(),
    subtitle="rat-snake: RNA, ADT and TCR single-cell analysis",
)

# 1 Summary -------------------------------------------------------------------
r.heading("Summary")
counts = status["status"].value_counts() if len(status) else pd.Series(dtype=int)
lines = [
    f"{len(status)} samples in {status['dataset'].nunique() if len(status) else 0} dataset(s): "
    + ", ".join(f"{n} {k}" for k, n in counts.items())
    + ".",
    f"Modalities: {', '.join(par.modalities)}. TCR clonotypes: {'yes' if par.tcr else 'no'}.",
]
if len(qc) and "cells_final" in qc.columns:
    lines.append(f"{int(qc['cells_final'].sum()):,} cells after preprocessing.")
lines.append(f"Populations: {', '.join(plan)}.")
r.bullets(lines)

# 2 Samples -------------------------------------------------------------------
r.heading("Samples")
r.text("Status of every sample in the sample tables. Excluded samples lack a required library; "
       "existing samples use a finished Cell Ranger run (`cellranger.existing_outs`).")
r.table(status.drop(columns=["cellranger_dir"], errors="ignore"), caption="Samples", max_rows=max_rows)

# 3 Cell Ranger ---------------------------------------------------------------
r.heading("Cell Ranger")
crm = read(inp.cellranger_metrics)
if crm.empty:
    r.text("_No Cell Ranger metrics._")
else:
    r.text("Cell Ranger multi metrics per sample (all metrics: `results/cellranger/metrics/`).")
    r.table(crm, caption="Cell Ranger metrics", max_rows=max_rows, max_cols=12)

# 4 Preprocessing -------------------------------------------------------------
r.heading("Preprocessing")
r.text("SoupX removes ambient RNA, DoubletFinder calls doublets, dsb normalises the ADT counts "
       "against empty droplets, and the RNA and ADT data are joined and filtered.")
if qc.empty:
    r.text("_No preprocessing summary._")
else:
    r.table(qc, caption="Cells and QC per sample", max_rows=max_rows)
for title, files, suffix in [
    ("Doublets", inp.doublet_qc, "_doublet_QC.png"),
    ("ADT background (dsb)", inp.adt_scatter, "_ADT_ScatterPlot.png"),
    ("RNA QC after filtering", inp.postfilt_vln, "_Postfilt_VlnPlot.png"),
]:
    r.heading(title, 2)
    for f in files:
        r.figure(f, caption=f"{title}: {sample_of(f, suffix)}")

# 5 TCR -----------------------------------------------------------------------
if par.tcr:
    r.heading("TCR clonotypes")
    r.text("Clonotypes from Cell Ranger's filtered contigs (scRepertoire combineTCR).")
    tcr_cols = [c for c in qc.columns if c in ("dataset", "sample", "tcr_cells", "clonotypes")
                or c in ("largest_clone",)]
    if tcr_cols:
        r.table(qc[tcr_cols], caption="Clonotypes per sample", max_rows=max_rows)
    for f in inp.clone_sizes:
        r.figure(f, caption=f"Clone sizes: {sample_of(f, '_clone_sizes.png')}")

# 6 Integrated object ---------------------------------------------------------
r.heading("Integrated object")
r.text("Samples are merged, integrated with Harmony (RNA and ADT) and joined with WNN.")
matches = read(inp.sample_matches)
if not matches.empty:
    r.table(matches, caption="Sample names matched to the metadata tables", max_rows=max_rows)
index = read(inp.metadata_plots)
for _, row in index.iterrows():
    r.figure(row["file"], caption=row["plot"].replace("__", ": ").replace("_", " "))

# 7 Populations ---------------------------------------------------------------
r.heading("Populations")
for pop, p in plan.items():
    r.heading(pop, 2)
    sizes = read(p["cluster_sizes"])
    if not sizes.empty:
        wide = sizes.assign(resolution=sizes["column"].str.replace(r".*_res\.", "", regex=True)).pivot_table(
            index="cluster", columns="resolution", values="cells", aggfunc="sum")
        r.table(wide.fillna(0).astype(int), caption=f"{pop}: cells per cluster at each resolution",
                max_rows=max_rows, index=True)
    r.figure(p["clustree"], caption=f"{pop}: clustree")
    for res in p["resolutions"]:
        key = str(res)
        r.heading(f"{pop}, resolution {res:g}", 3)
        for assay, path in p["markers"].get(key, {}).items():
            r.table(top_markers(path, top_n),
                    caption=f"{pop}, resolution {res:g}: top {top_n} {assay} markers per cluster "
                            f"(threshold set {opts.get('threshold_set', 'default')})",
                    max_rows=max_rows)
        for job in p["plots"].get(key, []):
            label = f"{pop}, resolution {res:g}, {job['assay']} {job['marker_set']}"
            r.figure(job["dotplot"], caption=f"{label}: dot plot")
            r.figure(job["markerplot"], caption=f"{label}: marker plot")

# 8 Settings ------------------------------------------------------------------
r.heading("Settings")
r.text("Main analysis options (the full config is in `config/config.yaml` of the output folder).")
r.text("```\n" + yaml.safe_dump(dict(par.settings), sort_keys=False, default_flow_style=None).strip() + "\n```")

prefix = snakemake.output.md[: -len(".md")]
formats = list(opts.get("formats", ["html", "docx", "notion"]))
written = r.render(prefix, formats=formats)
log_message(f"Report written: {', '.join(written.values())}", level="DONE")
