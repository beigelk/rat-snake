"""Final outputs of the annotate module for rule all."""

from rat_snake.context import context


def targets(config, paths):
    """Return the metadata table, each population's compiled object, clustree and every marker table.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    ctx = context(config, paths)
    pops = ctx.pops
    return [
        ctx.metadata_table,
        [pops.compiled_object(p) for p in pops.names],
        [pops.clustree(p) for p in pops.names],
        [pops.markers_filtered(p, r, a, t) for p, r, a, t in pops.filter_jobs()],
    ]
