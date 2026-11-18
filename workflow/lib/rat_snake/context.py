"""One object holding a run's samples and populations, shared by the Snakefile and every targets.py.

Reading the sample tables and checking finished Cell Ranger runs happens once
per config; ``context()`` hands back the same object to every caller.
"""

from rat_snake.populations import Populations
from rat_snake.samples import build_samples, cellranger_files, usable

_CACHE = {}


class RunContext:
    """Samples, populations and the shared paths of one run.

    Args:
        config: Snakemake config dict (already validated).
        paths: ``snake_core.paths.ProjectPaths`` for the run.
    """

    def __init__(self, config, paths):
        self.config = config
        self.P = paths
        self.samples = build_samples(config)
        self.usable = usable(self.samples)
        modalities = list(config.get("modalities", ["gex", "csp", "vdj"]))
        self.tcr = "vdj" in modalities and config.get("preprocess", {}).get("tcr", {}).get("enabled", True)
        self.integrated_object = paths.result("integrate", "harmony", "harmonized_seurat.qs2")
        self.metadata_table = paths.result("annotate", "metadata", "cell_metadata.qs2")
        self.pops = Populations(config, paths, self.integrated_object)

    def to_run(self):
        """``(dataset, sample)`` pairs Cell Ranger runs for."""
        return [(d, s) for d, s in self.usable if self.samples[d][s]["status"] == "run"]

    def cellranger_dir(self, dataset, sample):
        """Cell Ranger run folder: a finished run, or the one this pipeline makes."""
        info = self.samples[dataset][sample]
        if info["status"] == "existing":
            return info["run_dir"]
        return self.P.result("cellranger", "multi", dataset, sample)

    def cellranger_file(self, dataset, sample, key):
        """One Cell Ranger output of a sample (see ``samples.cellranger_files``)."""
        return cellranger_files(self.cellranger_dir(dataset, sample), sample)[key]


def context(config, paths):
    """Return the ``RunContext`` for this config, building it on first use."""
    key = id(config)
    if key not in _CACHE:
        _CACHE[key] = RunContext(config, paths)
    return _CACHE[key]
