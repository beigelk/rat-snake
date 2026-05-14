
rule cellranger_multi:
    input:
        RES_IDX['cellranger/configs']['outdir'] + "{dataset}/configfiles/{sample}_config.csv",
    output:
        touchfile = RES_IDX['cellranger/multi']['outdir'] + "touchfiles/{dataset}/{sample}_done.txt"
    params:
        outdir = RES_IDX['cellranger/multi']['outdir'] + "{dataset}",
        mem_gb = lambda wildcards, resources: resources.mem_mb // 1024
    log:
        LOG_IDX['cellranger/multi']['log'] + "{dataset}/{sample}.log"
    shell:
        """
        SCRATCH=/scr1/users/$USER/$SLURM_JOB_ID

        echo "=== Creating scratch dir ==="
        # Make scratch dir for sample
        mkdir -p $SCRATCH/{params.outdir}/{wildcards.sample}

        if [ $? -ne 0 ]; then
            echo "FAILED: could not create $SCRATCH"
            exit 1
        fi
        echo "OK: $SCRATCH created"

        echo "=== Running rule cellranger_multi ==="
        # Load CellRanger module
        module load CellRanger/9.0.1

        # Make final output dir (to copy to)
        mkdir -p {params.outdir}/{wildcards.sample}

        # Copy config to scratch
        cp --parents {input} $SCRATCH

        # Go to scratch dir where everything will be run (within dataset folder)
        cd $SCRATCH/{params.outdir}

        # Run cellranger multi
        cellranger multi \
            --id={wildcards.sample} \
            --csv=$SCRATCH/{input} \
            --output-dir=$SCRATCH/{params.outdir}/{wildcards.sample} \
            --localcores={threads} \
            --localmem={params.mem_gb} \
            --disable-ui \
            > $SLURM_SUBMIT_DIR/{log} 2>&1

        echo "=== Deleting SC_MULTI_CS from scratch ==="
        rm -rf $SCRATCH/{params.outdir}/{wildcards.sample}/SC_MULTI_CS
        if [ $? -ne 0 ]; then
            echo "FAILED: could not delete SC_MULTI_CS dir"
            exit 1
        fi
        echo "OK: deleted SC_MULTI_CS dir in $SCRATCH/{params.outdir}/{wildcards.sample}"

        echo "=== Copying results back to submit dir ==="
        cp -r $SCRATCH/{params.outdir}/{wildcards.sample} $SLURM_SUBMIT_DIR/{params.outdir}/{wildcards.sample}
        if [ $? -ne 0 ]; then
            echo "FAILED: could not copy back to $SLURM_SUBMIT_DIR"
            exit 1
        fi
        echo "OK: copied results to $SLURM_SUBMIT_DIR"

        touch {output.touchfile}

        echo "=== Cleaning up scratch ==="
        rm -rf $SCRATCH
        echo "OK: $SCRATCH removed"

        """

# PREVIOUS. ALL BUT SAMPLE 5 WERE RUN WITH THIS. 
# module load CellRanger/9.0.1

# rm -f {output.touchfile}
# rm -f {log}
# rm -f {params.outdir}/__{wildcards.sample}.mro
# rm -rf {params.outdir}/{wildcards.sample}

# mkdir -p {params.outdir}/{wildcards.sample}
# cd {params.outdir}

# cellranger multi \
#     --id={wildcards.sample} \
#     --csv={input} \
#     --output-dir={params.outdir}/{wildcards.sample} \
#     --localcores={threads} \
#     --localmem={params.mem_gb} \
#     --disable-ui \
#     > {log} 2>&1

# touch {output.touchfile}
# rm -rf {params.outdir}/{wildcards.sample}/SC_MULTI_CS


# ## TO ADD
# """
#     SCRATCH=/scr1/users/$USER/$SLURM_JOB_ID

#     echo "=== Creating scratch dir ==="
#     # Make scratch dir for sample
#     mkdir -p ${SCRATCH}/{params.outdir}/{wildcards.sample}

#     if [ $? -ne 0 ]; then
#         echo "FAILED: could not create ${SCRATCH}"
#         exit 1
#     fi
#     echo "OK: ${SCRATCH} created"

#     echo "=== Running rule cellranger_multi ==="
#     # Load CellRanger module
#     module load CellRanger/9.0.1

#     # Make final output dir (to copy to)
#     mkdir -p {params.outdir}/{wildcards.sample}

#     # Copy config to scratch
#     cp --parents {input} ${SCRATCH}

#     # Go to scratch dir where everything will be run (within dataset folder)
#     cd ${SCRATCH}/{params.outdir}

#     # Run cellranger multi
#     cellranger multi \
#         --id={wildcards.sample} \
#         --csv=${SCRATCH}/{input} \
#         --output-dir=${SCRATCH}/{params.outdir}/{wildcards.sample} \
#         --localcores={threads} \
#         --localmem={params.mem_gb} \
#         --disable-ui \
#         > $SLURM_SUBMIT_DIR/{log} 2>&1

#     echo "=== DELETING STORAGE-INTENSIVE SC_MULTI_CS DIR FROM SCRATCH ==="
#     rm -rf {params.outdir}/{wildcards.sample}/SC_MULTI_CS
#     if [ $? -ne 0 ]; then
#         echo "FAILED: could not deleted SC_MULTI_CS dir"
#         exit 1
#     fi
#     echo "OK: deleted SC_MULTI_CS dir in ${SCRATCH}/{params.outdir}/{wildcards.sample}"

#     echo "=== Copying results back to submit dir ==="
#     cp -r ${SCRATCH}/{params.outdir}/{wildcards.sample} $SLURM_SUBMIT_DIR/{params.outdir}/{wildcards.sample}
#     if [ $? -ne 0 ]; then
#         echo "FAILED: could not copy back to $SLURM_SUBMIT_DIR"
#         exit 1
#     fi
#     echo "OK: copied results to $SLURM_SUBMIT_DIR"

#     touch {output.touchfile}

#     echo "=== Cleaning up scratch ==="
#     rm -rf ${SCRATCH}
#     echo "OK: ${SCRATCH} removed"
# """

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
        