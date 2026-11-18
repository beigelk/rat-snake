"""Rule make_samplesheet: list one modality's FASTQs per library.

Moved from the rule's ``run:`` block unchanged in what it finds (the first
occurrence of the upper-case modality tag in a file name). An empty result now
fails with the folders that were searched instead of writing an empty table.

FASTQ names look like ``{anything}{MODALITY}[-_]{library}_L001_R1_001.fastq.gz``.
The library name is the part between the modality tag and the lane; the
full sample name (Cell Ranger's ``fastq_id``) is everything before the lane.

Output columns: SampleName, FullSampleName, InputDir, then {lane}_R1 and
{lane}_R2 for every lane seen.
"""

import csv
import glob
import os
import re
import sys
from collections import defaultdict

from snake_core.logging import log_message, redirect_to_log


def scan(dirs, filetype, modality):
    """Find one modality's FASTQs in the given folders.

    Args:
        dirs: Folders to search (every modality's input_dirs, since one folder
            often holds all libraries).
        filetype: File extension, e.g. ``fastq.gz``.
        modality: Modality tag as it appears in file names, e.g. ``GEX``.

    Returns:
        ``(data, full_names, sample_dirs)``: ``data[sample][lane][read] = path``,
        the full sample name and the folder of each library.
    """
    data = defaultdict(lambda: defaultdict(dict))
    full_names, sample_dirs = {}, {}
    tag = re.compile(rf"{re.escape(modality)}[-_]*")
    for d in dict.fromkeys(dirs):
        for f in sorted(glob.glob(os.path.join(d, f"*.{filetype}"))):
            name = os.path.basename(f)
            m = tag.search(name)
            if not m:
                continue
            after_mod = name[m.end():]
            lane_match = re.search(r"_L\d{3}_", after_mod)
            boundary = re.search(r"_L\d{3}_R[12]", name)
            read = re.search(r"_(R[12])_", after_mod)
            if not lane_match or not boundary or not read:
                continue
            sample = after_mod[: lane_match.start()].strip("-_")
            full_names.setdefault(sample, name[: boundary.start()].rstrip("-_"))
            sample_dirs.setdefault(sample, d)
            lane = lane_match.group(0).strip("_")
            data[sample][lane][read.group(1)] = f
    return data, full_names, sample_dirs


def write_samplesheet(path, data, full_names, sample_dirs):
    """Write the samplesheet TSV (one row per library)."""
    lanes = sorted({lane for lanes in data.values() for lane in lanes})
    header = ["SampleName", "FullSampleName", "InputDir"] + [f"{lane}_{r}" for lane in lanes for r in ("R1", "R2")]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        for sample in sorted(data):
            row = [sample, full_names[sample], sample_dirs[sample]]
            for lane in lanes:
                reads = data[sample].get(lane, {})
                row += [reads.get("R1", ""), reads.get("R2", "")]
            w.writerow(row)


def main(smk):
    """Snakemake entry point."""
    redirect_to_log(smk.log[0])
    ds = smk.config["datasets"][smk.wildcards.dataset]
    dirs = list(ds["input_dirs"].values())
    data, full_names, sample_dirs = scan(dirs, ds.get("filetype", "fastq.gz"), smk.wildcards.modality.upper())
    if not data:
        raise ValueError(f"No {smk.wildcards.modality.upper()} FASTQs (*.{ds.get('filetype', 'fastq.gz')}) found in {dirs}")
    write_samplesheet(smk.output.tsv, data, full_names, sample_dirs)
    log_message(f"{len(data)} {smk.wildcards.modality} libraries in {smk.output.tsv}", level="DONE")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
