"""Final outputs of the setup module for rule all."""

from rat_snake.context import context


def targets(config, paths):
    """Return the copied config, the sample status table and the samplesheets Cell Ranger needs.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    ctx = context(config, paths)
    datasets = sorted({d for d, _ in ctx.to_run()})
    sheets = [
        paths.result("setup", "samplesheet", f"samplesheet_{d}_{m}.tsv") for d in datasets for m in config["modalities"]
    ]
    return [paths.config_copy(), paths.result("setup", "sample_status.tsv"), sheets]
