"""Final outputs of the preprocess module for rule all."""

from rat_snake.context import context


def targets(config, paths):
    """Return each usable sample's Seurat object and TCR clones, and the QC summary table.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    ctx = context(config, paths)
    objs = [paths.result("preprocess", "construct", d, "objects", f"{s}_Seurat.qs2") for d, s in ctx.usable]
    clones = [paths.result("preprocess", "tcr_clones", d, "objects", f"{s}_combineTCR.qs2") for d, s in ctx.usable]
    return [objs, clones if ctx.tcr else [], paths.result("preprocess", "preprocess_summary.tsv")]
