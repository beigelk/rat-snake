rule soupx:
    input:
        filt_path = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix",
        raw_path  = RES_IDX['cellranger/multi']['outdir'] + "{dataset}/{sample}/outs/multi/count/raw_feature_bc_matrix"
    params:
        r_library   = config['resources']['r_library'],
        sample_name = "{sample}",
        npcs        = config['preprocess']['rna']['soupx']['npcs'],
        resolution  = config['preprocess']['rna']['soupx']['resolution']
    output:
        soupx_object  = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/objects/soupx__{sample}.qs2",
        soupx_figure  = RES_IDX['preprocess/rna/soupx']['outdir'] + "{dataset}/figures/soupx__{sample}.pdf"
    envmodules:
        "R/4.5.3"
    log: 
        logfile = LOG_IDX['preprocess/rna/soupx']['log'] + "{dataset}/{sample}.log"
    script:
        RULE_IDX['preprocess/rna/soupx']['r']

# rule run_soupx:
#     input:
#         script = 'Rscripts/run_soupx.R',
#         filt_path = lambda wildcards: f"{SAMPLE_DIR}/{wildcards.s}/outs/per_sample_outs/{wildcards.s}/count/sample_filtered_feature_bc_matrix",
#         raw_path = lambda wildcards: f"{SAMPLE_DIR}/{wildcards.s}/outs/multi/count/raw_feature_bc_matrix"
#     params:
#         sample_name = '{s}'
#     output:
#         soupx_object = outdir + 'output/soupx/objects/{s}_SoupX.qs2',
#         soupx_figures = outdir + 'output/soupx/figures/{s}_SoupX.pdf'
#     log: 
#         logfile = logdir + 'run_soupx/run_soupx_{s}.log'
#     resources:
#         tmpdir = 'tmp/'
#     shell:
#         "( \
#         Rscript \
#         {input.script} \
#         {input.filt_path} \
#         {input.raw_path} \
#         {params.sample_name} \
#         {output.soupx_object} \
#         {output.soupx_figures} \
#         ) \
#         > {log.logfile} 2>&1"