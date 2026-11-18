# ──────────────────────────────────────────────
# downstream_step.smk
# ──────────────────────────────────────────────

wildcard_constraints:
    compile_target = "|".join(COMPILE_TARGETS),

def downstream_step_input(wc):
    """Pulls the compiled object for this target (main, epithelial, lymphocytes, ...)."""
    return RES_IDX['annotate/compile']['outdir'] + "objects/" + "compiled__" + wc.compile_target + ".qs2"

rule downstream_step:
    input:
        compiled_obj = downstream_step_input
    params:
        compile_target = lambda wc: wc.compile_target,
        r_library      = config['resources']['r_library']
    output:
        result_obj = RES_IDX['annotate/downstream_step']['outdir'] + "objects/" +
            "downstream_step__{compile_target}.qs2"
    envmodules:
        "R/4.5.3"
    log:
        logfile = LOG_IDX['annotate/downstream_step']['log'] + "downstream_step__{compile_target}.log"
    script:
        RULE_IDX['annotate/downstream_step']['r']

def complete_downstream_step_outputs(wildcards):
    targets = []
    for target in COMPILE_TARGETS:
        targets.append(
            RES_IDX['annotate/downstream_step']['outdir'] + "objects/" + "downstream_step__" + target + ".qs2"
        )
    return targets