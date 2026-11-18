"""Samples: who they are, which libraries they have and where their Cell Ranger output is.

Each dataset in the config has a ``sample_ids`` table (tab-separated) with one
row per sample and one column per modality (``gex``, ``csp``, ``vdj``) holding
that library's FASTQ sample name, for example::

    gex                 csp                 vdj
    25_13_1205_S25      25_13_1205_S41      25_13_1205_S57

A sample is named by the name parts its libraries share (``25_13_1205``), or
by an optional ``sample`` column. Each sample then gets one status:

``run``       Cell Ranger multi runs from its FASTQs.
``existing``  ``cellranger.existing_outs`` points at a finished Cell Ranger
              run, which is used as it is; Cell Ranger never runs for it.
``excluded``  A required library is missing from the table; the sample is
              left out (and listed in ``results/setup/sample_status.tsv``).
"""

import csv
import glob
import os
import re

# Cell Ranger multi feature type of each modality
FEATURE_TYPES = {"gex": "Gene Expression", "csp": "Antibody Capture", "vdj": "VDJ-T"}

# Sample-number tokens added by the sequencer, e.g. "S25"
_SAMPLE_NUMBER = re.compile(r"^S\d+$")


def shared_tokens(names):
    """Name a sample by the "_"-separated parts all its library names share.

    Sequencer sample numbers (``S`` followed by digits) are dropped. The old
    Snakefile dropped every part starting with "S", which also removed parts
    such as "Spleen".

    Args:
        names: Library names of one sample, e.g. ``["25_13_1205_S25", "25_13_1205_S41"]``.

    Returns:
        The shared parts joined by "_", e.g. ``"25_13_1205"``.
    """
    split = [n.split("_") for n in names]
    if not split:
        return ""
    shared = [t for t in split[0] if all(t in parts for parts in split[1:]) and not _SAMPLE_NUMBER.match(t)]
    return "_".join(shared)


def read_sample_table(path, modalities=("gex", "csp", "vdj")):
    """Read one dataset's sample table.

    Args:
        path: Tab-separated file with a column per modality and an optional
            ``sample`` column.
        modalities: Modality columns to read.

    Returns:
        Dict ``{sample: {modality: library name}}``, in file order. Empty
        cells are left out, so a sample may lack a modality.

    Raises:
        ValueError: A column is missing, a sample name is empty or two rows
            get the same name.
    """
    samples = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        missing = [m for m in modalities if m not in cols]
        if missing:
            raise ValueError(f"{path}: missing column(s) {missing}; found {cols}")
        for i, row in enumerate(reader, start=2):
            libs = {m: (row.get(m) or "").strip() for m in modalities}
            libs = {m: v for m, v in libs.items() if v}
            if not libs:
                continue
            name = (row.get("sample") or "").strip() or shared_tokens(list(libs.values()))
            if not name:
                raise ValueError(f"{path}, line {i}: the libraries share no name part; add a 'sample' column")
            if name in samples:
                raise ValueError(f"{path}, line {i}: sample name '{name}' is used twice; add a 'sample' column")
            samples[name] = libs
    return samples


def resolve_run_dir(path):
    """Return the Cell Ranger run folder (the one holding ``outs/``) for a path.

    Accepts the run folder itself or its ``outs/`` folder.

    Raises:
        FileNotFoundError: Neither form exists.
    """
    path = os.path.abspath(os.path.expanduser(str(path)))
    if os.path.isdir(os.path.join(path, "outs")):
        return path
    if os.path.basename(path.rstrip("/")) == "outs" and os.path.isdir(path):
        return os.path.dirname(path.rstrip("/"))
    raise FileNotFoundError(f"No Cell Ranger outs/ folder at {path}")


def existing_run_dir(dataset, sample, cellranger_cfg):
    """Find a finished Cell Ranger run for a sample, if the config names one.

    ``cellranger.existing_outs.samples`` maps sample names (or
    ``dataset/sample``) to run folders. Otherwise
    ``cellranger.existing_outs.root`` is searched for ``{root}/{dataset}/{sample}``
    and ``{root}/{sample}``.

    Returns:
        The run folder, or None when the sample should be run.
    """
    existing = (cellranger_cfg or {}).get("existing_outs") or {}
    by_sample = existing.get("samples") or {}
    for key in (f"{dataset}/{sample}", sample):
        if key in by_sample:
            return resolve_run_dir(by_sample[key])
    root = existing.get("root")
    if root:
        for cand in (os.path.join(root, dataset, sample), os.path.join(root, sample)):
            if os.path.isdir(os.path.join(cand, "outs")):
                return os.path.abspath(cand)
    return None


def cellranger_files(run_dir, sample):
    """Paths of the Cell Ranger multi outputs the pipeline reads.

    ``per_sample_outs/{sample}`` is used when it exists (runs made by this
    pipeline use ``--id {sample}``); otherwise the only folder under
    ``per_sample_outs`` (finished runs made elsewhere).

    Args:
        run_dir: Cell Ranger run folder (holds ``outs/``).
        sample: Sample name.

    Returns:
        Dict with ``raw_matrix``, ``filtered_matrix``, ``contigs``,
        ``metrics`` and ``web_summary``.
    """
    per_sample_root = os.path.join(run_dir, "outs", "per_sample_outs")
    per_sample = os.path.join(per_sample_root, sample)
    if not os.path.isdir(per_sample) and os.path.isdir(per_sample_root):
        subdirs = sorted(d for d in glob.glob(os.path.join(per_sample_root, "*")) if os.path.isdir(d))
        if len(subdirs) == 1:
            per_sample = subdirs[0]
    return {
        "raw_matrix": os.path.join(run_dir, "outs", "multi", "count", "raw_feature_bc_matrix"),
        "filtered_matrix": os.path.join(per_sample, "count", "sample_filtered_feature_bc_matrix"),
        "contigs": os.path.join(per_sample, "vdj_t", "filtered_contig_annotations.csv"),
        "metrics": os.path.join(per_sample, "metrics_summary.csv"),
        "web_summary": os.path.join(per_sample, "web_summary.html"),
    }


def check_existing_run(run_dir, sample, need_vdj):
    """List the outputs a finished Cell Ranger run is missing.

    Returns:
        A list of missing paths (empty when the run is usable).
    """
    files = cellranger_files(run_dir, sample)
    needed = ["raw_matrix", "filtered_matrix", "metrics"] + (["contigs"] if need_vdj else [])
    return [files[k] for k in needed if not os.path.exists(files[k])]


def build_samples(config):
    """Read every dataset's samples and give each a status.

    Args:
        config: Snakemake config dict. Reads ``datasets``, ``modalities``,
            ``cellranger.existing_outs`` and ``preprocess.tcr.enabled``.

    Returns:
        Dict ``{dataset: {sample: info}}``. ``info`` has ``libraries``
        (modality -> library name), ``status`` (run / existing / excluded),
        ``reason`` and, for usable samples, ``run_dir`` for existing runs.

    Raises:
        ValueError: A finished run named in the config lacks required outputs.
    """
    modalities = list(config.get("modalities", ["gex", "csp", "vdj"]))
    need_vdj = "vdj" in modalities and config.get("preprocess", {}).get("tcr", {}).get("enabled", True)
    out = {}
    for dataset, dcfg in config["datasets"].items():
        if not dcfg.get("enabled", True):
            continue
        out[dataset] = {}
        if not os.path.isfile(dcfg["sample_ids"]):
            raise FileNotFoundError(f"datasets.{dataset}.sample_ids: no file at {dcfg['sample_ids']}")
        for sample, libs in read_sample_table(dcfg["sample_ids"], modalities).items():
            info = {"libraries": libs, "status": "run", "reason": "", "run_dir": None}
            run_dir = existing_run_dir(dataset, sample, config.get("cellranger"))
            if run_dir:
                missing = check_existing_run(run_dir, sample, need_vdj)
                if missing:
                    raise ValueError(
                        f"cellranger.existing_outs for {dataset}/{sample} ({run_dir}) is missing: " + ", ".join(missing)
                    )
                info.update(status="existing", reason=f"Cell Ranger outputs from {run_dir}", run_dir=run_dir)
            else:
                lacking = [m for m in modalities if m not in libs]
                if lacking:
                    info.update(status="excluded", reason=f"no {', '.join(lacking)} library in the sample table")
            out[dataset][sample] = info
    return out


def usable(samples):
    """``(dataset, sample)`` pairs that go on to preprocessing, in table order."""
    return [(d, s) for d, ss in samples.items() for s, i in ss.items() if i["status"] != "excluded"]
