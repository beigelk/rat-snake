# ==============================================================================
# MODULE: setup
# Record what a run used (the config and every sample's status) and list each
# dataset's FASTQs per modality for Cell Ranger.
# ==============================================================================

SETUP_SCRIPTS = os.path.join(MODULES_DIR, "setup", "scripts")

localrules:
    copy_config,
    sample_status,


rule copy_config:
    """Copy the run's config (after defaults are filled in) into the output folder."""
    output:
        config=P.config_copy(),
    log:
        P.log("setup", "copy_config"),
    run:
        import yaml

        with open(output.config, "w") as fh:
            yaml.safe_dump(config, fh, sort_keys=False)
        log_message(f"Copied config to {output.config}", level="INFO")


rule sample_status:
    """Write every sample's libraries and status (run, existing or excluded, and why)."""
    output:
        tsv=P.result("setup", "sample_status.tsv"),
    log:
        P.log("setup", "sample_status"),
    run:
        import csv

        with open(output.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["dataset", "sample", *MODALITIES, "status", "reason", "cellranger_dir"])
            for d, ss in SAMPLES.items():
                for s, info in ss.items():
                    crdir = CTX.cellranger_dir(d, s) if info["status"] != "excluded" else ""
                    libs = [info["libraries"].get(m, "") for m in MODALITIES]
                    w.writerow([d, s, *libs, info["status"], info["reason"], crdir])


rule make_samplesheet:
    """List one dataset's FASTQs for one modality, one row per library with its lanes and reads.

    Searches every folder in `datasets.<dataset>.input_dirs` for files with
    the upper-case modality tag (GEX, CSP, VDJ) in their name.
    """
    output:
        tsv=P.result("setup", "samplesheet", "samplesheet_{dataset}_{modality}.tsv"),
    log:
        P.log("setup", "make_samplesheet", "{dataset}", "{modality}"),
    script:
        os.path.join(SETUP_SCRIPTS, "make_samplesheet.py")
