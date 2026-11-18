"""Final outputs of the report module for rule all."""

FILES = {"html": "report.html", "docx": "report.docx", "notion": "report_notion.docx"}


def targets(config, paths):
    """Return the report in every format of ``report.formats``.

    Args:
        config: Snakemake config dict.
        paths: snake_core.paths.ProjectPaths for this run.
    """
    return [paths.report(FILES[f]) for f in config["report"]["formats"]]
