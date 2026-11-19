# ==============================================================================
# MODULE: cellranger
# Cell Ranger multi per sample (GEX + CSP + VDJ), or a finished run from
# `cellranger.existing_outs`, then one metrics table and a zip of the web
# summaries.
# ==============================================================================

CR_SCRIPTS = os.path.join(MODULES_DIR, "cellranger", "scripts")
CR_CFG = config["cellranger"]
CR_RUN_DIR = P.result("cellranger", "multi", "{dataset}", "{sample}")

localrules:
    make_cellranger_config,
    collect_web_summaries,


def samplesheets_for(wc):
    """The dataset's samplesheet of every modality, by modality."""
    return {m: P.result("setup", "samplesheet", f"samplesheet_{wc.dataset}_{m}.tsv") for m in MODALITIES}


def cellranger_runs(wc=None):
    """Cell Ranger run folder of every usable sample (run here or finished earlier)."""
    return [CTX.cellranger_dir(d, s) for d, s in USABLE]


rule make_cellranger_config:
    """Write one sample's Cell Ranger multi config from the samplesheets.

    Libraries come from the dataset's `sample_ids` table; their FASTQ ids and
    folders from the samplesheets. Fails if a library has no FASTQs.
    """
    input:
        unpack(samplesheets_for),
    output:
        csv=P.result("cellranger", "configs", "{dataset}", "{sample}_config.csv"),
    params:
        modalities=MODALITIES,
        libraries=lambda wc: SAMPLES[wc.dataset][wc.sample]["libraries"],
    log:
        P.log("cellranger", "make_cellranger_config", "{dataset}", "{sample}"),
    script:
        os.path.join(CR_SCRIPTS, "make_cellranger_config.py")


rule cellranger_multi:
    """Run Cell Ranger multi for one sample.

    Runs in `cellranger.scratch_dir` when set (copying the run back), else
    next to the output. `cellranger.path` or `software.envmodules.cellranger`
    provides Cell Ranger. Samples with a finished run in
    `cellranger.existing_outs` never reach this rule.
    """
    input:
        csv=rules.make_cellranger_config.output.csv,
    output:
        run_dir=directory(CR_RUN_DIR),
    params:
        script=os.path.join(CR_SCRIPTS, "run_cellranger_multi.sh"),
        scratch=CR_CFG.get("scratch_dir") or "",
        keep_cs=str(CR_CFG.get("keep_sc_multi_cs", False)).lower(),
        cellranger=CR_CFG.get("path") or "",
        mem_gb=CR_CFG.get("mem_gb", 64),
    threads: CR_CFG.get("threads", 16)
    resources:
        mem_mb=CR_CFG.get("mem_gb", 64) * 1024 + 4096,
        # Slurm: ask for as many CPUs as Cell Ranger uses. Without this, the
        # profile's default cpus_per_task (1) wins over threads and all
        # --localcores share one CPU.
        cpus_per_task=CR_CFG.get("threads", 16),
        cellranger_jobs=1,
    envmodules:
        *tool_envmodules("cellranger", config),
    log:
        P.log("cellranger", "cellranger_multi", "{dataset}", "{sample}"),
    benchmark:
        P.benchmark("cellranger", "cellranger_multi", "{dataset}", "{sample}")
    shell:
        "bash {params.script} {input.csv:q} {wildcards.sample:q} {output.run_dir:q} {threads} {params.mem_gb}"
        " {params.scratch:q} {params.keep_cs} {params.cellranger:q} > {log:q} 2>&1"


rule cellranger_metrics:
    """Combine every sample's Cell Ranger metrics_summary.csv into one table (sample x metric)."""
    input:
        runs=cellranger_runs,
    output:
        tsv=P.result("cellranger", "metrics", "cellranger_metrics.tsv"),
        long=P.result("cellranger", "metrics", "cellranger_metrics_long.tsv"),
    params:
        metrics={s: CTX.cellranger_file(d, s, "metrics") for d, s in USABLE},
    log:
        P.log("cellranger", "cellranger_metrics"),
    conda:
        r_env("cellranger").conda
    container:
        r_env("cellranger").container
    script:
        os.path.join(CR_SCRIPTS, "cellranger_metrics.R")


rule collect_web_summaries:
    """Zip every sample's Cell Ranger web_summary.html, renamed to {dataset}__{sample}_web_summary.html."""
    input:
        runs=cellranger_runs,
    output:
        zip=P.result("cellranger", "web_summaries", f"{PROJECT}_CellRanger_web_summaries.zip"),
    params:
        files={f"{d}__{s}_web_summary.html": CTX.cellranger_file(d, s, "web_summary") for d, s in USABLE},
    log:
        P.log("cellranger", "collect_web_summaries"),
    run:
        import zipfile

        with zipfile.ZipFile(output.zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for name, path in params.files.items():
                if os.path.exists(path):
                    zf.write(path, arcname=name)
                else:
                    log_message(f"No web summary at {path}", level="WARN")
