
rule run_cellranger_multi:
    input:
        RES_IDX['cellranger/configs']['outdir'] + "{dataset}/configfiles/{sample}_config.csv",
    output:
        RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv"
    params:
        outdir=RES_IDX['cellranger/multi']['outdir'] + "{dataset}",
        mem_gb=lambda wildcards, resources: resources.mem_mb // 1024
    threads: 16
    resources:
        mem_mb=64000,
        runtime=720
    log:
        LOG_IDX['cellranger/multi']['log'] + "{dataset}/{sample}.log"
    shell:
        """
        module load CellRanger/9.0.1

        rm -rf {params.outdir}/{wildcards.sample}

        cellranger multi \
            --id={wildcards.sample} \
            --csv={input} \
            --output-dir={params.outdir}/{wildcards.sample} \
            --localcores={threads} \
            --localmem={params.mem_gb} \
            --disable-ui \
            > {log} 2>&1
        """



# def complete_cellranger_inputs(wildcards):
#     # forces snakemake to wait until all make_configs jobs are done
#     configs_dir = RES_IDX['cellranger/configs']['outdir'] + f"{wildcards.dataset}/"
#     # trigger checkpoint evaluation
#     checkpoint_output = checkpoints.make_configs.get(**wildcards).output.skipped

#     complete = []
#     for dataset, samples in SAMPLES.items():
#         for sample_key in samples:
#             skipped_path = (
#                 RES_IDX['cellranger/configs']['outdir'] + f"{dataset}/tracking/{sample_key}_record.tsv"
#             )
#             with open(skipped_path) as f:
#                 rows = list(csv.DictReader(f, delimiter="\t"))
#             if not rows:
#                 complete.append((dataset, sample_key))

#     return expand(
#         RES_IDX['cellranger/configs']['outdir'] + "{dataset}/configfiles/{sample}_config.csv",
#         zip,
#         dataset=[d for d, s in complete],
#         sample=[s for d, s in complete]
#     )


# rule run_cellranger_multi:
#     input:
#         complete_cellranger_inputs
#         # RES_IDX['cellranger/configs']['outdir'] + "cellranger_configs/{dataset}/{sample}_config.csv",
#     output:
#         RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv"
#     threads: 16
#     resources:
#         mem_mb=64000
#     shell:
#         "cellranger multi"
#         " --id={wildcards.sample}"
#         " --csv={input}"
#         " --output-dir=" + RES_IDX['cellranger/multi']['outdir'] + "{wildcards.dataset}/{wildcards.sample}"
#         " --localcores={threads}"
#         " --localmem={resources.mem_mb}"

# # rule run_cellranger_multi:
# #     input:
# #         complete_cellranger_inputs
# #         # config=RES_IDX['cellranger/configs']['outdir'] + "{dataset}/configfiles/{sample}_config.csv",
# #         # skipped=RES_IDX['cellranger/configs']['outdir'] + "{dataset}/tracking/{sample}_record.tsv",
# #     output:
# #         RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv"
# #     run:
# #         import csv

# #         # check skipped file
# #         with open(input.skipped) as f:
# #             reader = csv.DictReader(f, delimiter="\t")
# #             skipped_rows = list(reader)

# #         if skipped_rows:
# #             print(f"Skipping cellranger multi for {wildcards.sample} — missing modalities:")
# #             for row in skipped_rows:
# #                 print(f"  {row['modality']}: {row['reason']}")
# #         else:
# #             shell(
# #                 "cellranger multi"
# #                 " --id={wildcards.sample}"
# #                 " --csv={input.config}"
# #                 " --output-dir=" + RES_IDX['cellranger/multi']['outdir'] + "{wildcards.dataset}/{wildcards.sample}"
# #                 " --localcores={threads}"
# #                 " --localmem={resources.mem_mb}"
# #             )




# rule cellranger_multi:
#     input:
#     output:
#     shell:
#         """
#         cellranger multi --localmem=200 --id=F422-3_S10 --csv=config_files/F422-3_S10_config.csv
#         """
        