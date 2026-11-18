"""Final outputs of the integrate module for rule all."""

from rat_snake.context import context


def targets(config, paths):
    """Return the integrated (Harmony + WNN) object.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    return [context(config, paths).integrated_object]
