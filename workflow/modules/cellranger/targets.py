"""Final outputs of the cellranger module for rule all."""


def targets(config, paths):
    """Return the metrics table and the web summaries zip (which need every sample's run).

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    return [
        paths.result("cellranger", "metrics", "cellranger_metrics.tsv"),
        paths.result("cellranger", "web_summaries", f"{paths.project}_CellRanger_web_summaries.zip"),
    ]
