rule collect_web_summaries:
    input:
        complete_cellranger_outputs
    params:
        web_summaries = expand(
            RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
            zip,
            dataset=[d for d, s in get_complete_samples()],
            sample=[s for d, s in get_complete_samples()]
        )
    output:
        zip = RES_IDX['utils/collect_web_summaries']['outdir'] + PROJECT + "_CellRanger_web_summaries.zip"
    log:
        LOG_IDX['utils/collect_web_summaries']['log'] + PROJECT + "_zip.log"
    shell:
        """
        zip {output.zip} {params.web_summaries} > {log} 2>&1
        """