"""Final outputs of the visualization module for rule all."""

from rat_snake.context import context


def targets(config, paths):
    """Return the metadata plot index and every dot plot and marker plot.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    pops = context(config, paths).pops
    return [
        paths.result("visualization", "metadata", "metadata_plots.tsv"),
        [pops.dotplot(*job) for job in pops.plot_jobs()],
        [pops.markerplot(*job) for job in pops.plot_jobs()],
    ]
