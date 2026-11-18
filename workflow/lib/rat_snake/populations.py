"""Populations: the clustering, subsetting and marker plan in ``annotate.populations``.

A population is a set of cells clustered at several resolutions. A *root*
population (no ``subset`` block) is the whole integrated object. Any other
population is cut from a parent population's clusters at one resolution and
re-clustered (``subset: {base_name, base_resolution, base_clusters}``), so
populations form a tree, e.g. main -> lymphocytes -> CD8s.

Every file name used by the annotate and visualization modules is built here,
so the rules, the target lists and the report agree on them. The path methods
also take wildcard placeholders (``"{population}"``, ``"{res}"``) to build rule
output patterns. File names leave out the graph name; the cluster column
names inside the objects keep it (``cluster__{pop}__{graph}_res.{res}``). Resolutions are
written with ``format_res`` (``0.0 -> "0"``, ``0.10 -> "0.1"``), which matches
how R prints the same number in the cluster column names.
"""

import os

import numpy as np


def format_res(value):
    """Write a resolution the way R's ``paste0()`` does.

    ``0.0 -> "0"``, ``1.0 -> "1"``, ``0.1 -> "0.1"``, ``0.07 -> "0.07"``.
    A wildcard placeholder such as ``"{res}"`` is returned unchanged, so the
    path methods below also build rule output patterns.
    """
    if isinstance(value, str) and value.startswith("{"):
        return value
    value = round(float(value), 4)
    if value == int(value):
        return str(int(value))
    return f"{value:g}"


def parse_clusters(clusters):
    """Read ``base_clusters``: ``"0, 1, 8"``, ``4``, ``[0, 1]`` or ``"all"``.

    Returns:
        A list of ints, or the string ``"all"`` (keep every cluster).
    """
    if isinstance(clusters, str) and clusters.strip().lower() == "all":
        return "all"
    if isinstance(clusters, str):
        return [int(c.strip()) for c in clusters.split(",") if c.strip()]
    if isinstance(clusters, (int, float)):
        return [int(clusters)]
    return [int(c) for c in clusters]


def expand_resolutions(spec):
    """Resolutions from one scheme: ``{from, to, by}`` or ``{values: [...]}``.

    ``to`` is included. Values are rounded to 4 decimals.
    """
    if "values" in spec:
        return [round(float(v), 4) for v in spec["values"]]
    missing = {"from", "to", "by"} - set(spec)
    if missing:
        raise KeyError(f"resolution scheme needs from, to and by (missing {sorted(missing)}): {spec}")
    vals = np.arange(spec["from"], spec["to"] + spec["by"] / 2, spec["by"])
    return [round(float(v), 4) for v in vals]


class Populations:
    """The populations of one run and the paths of everything made for them.

    Args:
        config: Snakemake config dict (reads ``annotate.populations``).
        paths: ``snake_core.paths.ProjectPaths`` for the run.
        root_object: Path of the object root populations are clustered on
            (the integrated object).
    """

    def __init__(self, config, paths, root_object):
        self.P = paths
        self.project = paths.project
        self.root_object = root_object
        pops = config.get("annotate", {}).get("populations") or []
        self.pops = {}
        for p in pops:
            name = str(p["name"]).strip()
            if name in self.pops:
                raise ValueError(f"annotate.populations: '{name}' is listed twice")
            self.pops[name] = {**p, "name": name}
        for name, p in self.pops.items():
            if "subset" in p:
                base = p["subset"]["base_name"]
                if base not in self.pops:
                    raise ValueError(f"Population '{name}' subsets unknown population '{base}'")
                res = round(float(p["subset"]["base_resolution"]), 4)
                if res not in self.resolutions(base):
                    raise ValueError(
                        f"Population '{name}' uses resolution {format_res(res)} of '{base}', "
                        f"which '{base}' is not clustered at ({', '.join(map(format_res, self.resolutions(base)))})"
                    )

    # -- the plan --------------------------------------------------------------
    @property
    def names(self):
        """Population names in config order."""
        return list(self.pops)

    def is_root(self, pop):
        """True when the population is the whole integrated object."""
        return "subset" not in self.pops[pop]

    def graph(self, pop):
        """Graph the population is clustered on (default ``wsnn_harmony``)."""
        return self.pops[pop].get("graph_name", "wsnn_harmony")

    def resolutions(self, pop):
        """Sorted resolutions of a population, from all its schemes."""
        res = set()
        for spec in (self.pops[pop].get("resolutions") or {}).values():
            res.update(expand_resolutions(spec))
        return sorted(res)

    def subset_spec(self, pop):
        """The parent, its resolution and the clusters a population keeps."""
        s = self.pops[pop]["subset"]
        return {
            "base_name": s["base_name"],
            "base_graph_name": self.graph(s["base_name"]),
            "base_resolution": format_res(s["base_resolution"]),
            "base_clusters": parse_clusters(s["base_clusters"]),
        }

    def wnn(self, pop, defaults):
        """WNN settings used when a subset is re-analysed (population ``wnn`` over ``defaults``)."""
        return {**(defaults or {}), **(self.pops[pop].get("wnn") or {})}

    def assays(self, pop):
        """Assays markers are found for (``findallmarkers.assays``)."""
        return list((self.pops[pop].get("findallmarkers") or {}).get("assays", []))

    def threshold_sets(self, pop, assay):
        """``{name: {lfc_threshold, padj_threshold, pct1_threshold}}`` for an assay."""
        return dict(((self.pops[pop].get("findallmarkers") or {}).get("threshold_sets") or {}).get(assay, {}))

    def marker_sets(self, pop, assay):
        """``{name: spec}`` plot marker sets for an assay (``plots.marker_sets``)."""
        return dict(((self.pops[pop].get("plots") or {}).get("marker_sets") or {}).get(assay, {}))

    # -- paths -----------------------------------------------------------------
    def _res(self, module, pop, *parts):
        return self.P.result(module, "populations", pop, *parts)

    def cluster_column(self, pop, res):
        """Metadata column holding a population's clusters at one resolution."""
        return f"cluster__{pop}__{self.graph(pop)}_res.{format_res(res)}"

    def cluster_tsv(self, pop, res):
        """Cluster table (barcode -> cluster) of one resolution."""
        return self._res("annotate", pop, "cluster", f"cluster__{pop}_res.{format_res(res)}.tsv")

    def subset_object(self, pop):
        """Re-analysed object of a subset population."""
        return self._res("annotate", pop, "subset", f"subset__{pop}.qs2")

    def input_object(self, pop):
        """Object a population is clustered on."""
        return self.root_object if self.is_root(pop) else self.subset_object(pop)

    def compiled_object(self, pop):
        """Population object with every resolution's clusters and the sample metadata."""
        return self._res("annotate", pop, "compile", f"{pop}_seurat.qs2")

    def cluster_sizes(self, pop):
        """Cells per cluster at every resolution (written by compile)."""
        return self._res("annotate", pop, "compile", f"{pop}_cluster_sizes.tsv")

    def markers_raw(self, pop, res, assay):
        """FindAllMarkers table of one resolution and assay."""
        return self._res(
            "annotate", pop, "findallmarkers", f"{self.project}__{pop}__{assay}_res.{format_res(res)}__markers_raw.tsv"
        )

    def markers_filtered(self, pop, res, assay, tset):
        """Markers passing one threshold set."""
        return self._res(
            "annotate",
            pop,
            "filtermarkers",
            f"{self.project}__{pop}__{assay}_res.{format_res(res)}__{tset}__markers_filtered.tsv",
        )

    def clustree(self, pop):
        """clustree figure (PNG; a PDF sits next to it)."""
        return self._res("annotate", pop, "clustree", f"{self.project}__{pop}__clustree.png")

    def plot_prefix(self, pop, res, assay, mset):
        """Shared start of a population's plot file names."""
        return f"{self.project}__{pop}__{assay}_res.{format_res(res)}__{mset}"

    def dotplot(self, pop, res, assay, mset):
        """Dot plot (PNG) of one marker set."""
        return self._res("visualization", pop, "dotplot", f"{self.plot_prefix(pop, res, assay, mset)}__DotPlot.png")

    def markerplot(self, pop, res, assay, mset):
        """Combined UMAP, feature and violin plot (PNG) of one marker set."""
        return self._res(
            "visualization", pop, "markerplots", f"{self.plot_prefix(pop, res, assay, mset)}__ComboPlot.png"
        )

    # -- targets ---------------------------------------------------------------
    def marker_jobs(self):
        """``(pop, res, assay)`` for every FindAllMarkers job."""
        return [(p, r, a) for p in self.names for r in self.resolutions(p) for a in self.assays(p)]

    def filter_jobs(self):
        """``(pop, res, assay, threshold_set)`` for every marker filter job."""
        return [(p, r, a, t) for p, r, a in self.marker_jobs() for t in self.threshold_sets(p, a)]

    def plot_jobs(self):
        """``(pop, res, assay, marker_set)`` for every plot job."""
        return [(p, r, a, m) for p, r, a in self.marker_jobs() for m in self.marker_sets(p, a)]

    def plot_features_input(self, pop, res, assay, mset):
        """File the features of a marker set come from.

        ``topn`` sets read a filtered marker table; ``genelist`` sets read the
        listed file (one feature per line).
        """
        spec = self.marker_sets(pop, assay)[mset]
        if spec["type"] == "topn":
            tset = spec.get("threshold_set", "default")
            if tset not in self.threshold_sets(pop, assay):
                raise ValueError(f"Marker set {pop}/{assay}/{mset} uses unknown threshold set '{tset}'")
            return self.markers_filtered(pop, res, assay, tset)
        if spec["type"] == "genelist":
            return os.path.abspath(spec["file"])
        raise ValueError(f"Marker set {pop}/{assay}/{mset}: type must be topn or genelist (got {spec['type']!r})")
