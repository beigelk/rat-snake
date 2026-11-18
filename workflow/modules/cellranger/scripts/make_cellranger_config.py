"""Rule make_cellranger_config: write one sample's Cell Ranger multi config CSV.

From the old ``make_configs`` checkpoint. Samples missing a library are now
left out before the run (see rat_snake.samples), so this rule only runs for
complete samples and fails loudly if a library has no FASTQs, instead of
writing a config that skips it. Sections are written only for the
modalities in ``modalities``.
"""

import csv
import os
import re

from snake_core.logging import log_message, redirect_to_log

# Cell Ranger multi feature type of each modality
FEATURE_TYPES = {"gex": "Gene Expression", "csp": "Antibody Capture", "vdj": "VDJ-T"}


def load_samplesheet(path):
    """Read a samplesheet: library name -> Cell Ranger ``fastq_id`` and FASTQ folder.

    ``fastq_id`` is the full sample name without the sequencer sample number
    (``_S25``), as Cell Ranger expects.
    """
    index = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            index[row["SampleName"]] = {
                "fastq_id": re.sub(r"_S\d+$", "", row["FullSampleName"]),
                # Absolute: Cell Ranger runs in its own (scratch) folder
                "fastqs": os.path.abspath(row["InputDir"]),
            }
    return index


def config_lines(cr, libraries, modalities):
    """Lines of a Cell Ranger multi config.

    Args:
        cr: The config's ``cellranger`` block.
        libraries: ``[{fastq_id, fastqs, feature_types}]``.
        modalities: Modalities in this run.
    """
    lines = [
        "[gene-expression]",
        f"reference,{os.path.abspath(cr['gex_reference'])}",
        f"create-bam,{str(cr.get('create_bam', False)).lower()}",
        f"check-library-compatibility,{str(cr.get('check_library_compatibility', False)).lower()}",
        "",
    ]
    if "csp" in modalities:
        lines += ["[feature]", f"reference,{os.path.abspath(cr['csp_reference'])}", ""]
    if "vdj" in modalities:
        lines += ["[vdj]", f"reference,{os.path.abspath(cr['vdj_reference'])}", ""]
    lines += ["[libraries]", "fastq_id,fastqs,lanes,feature_types"]
    lanes = cr.get("lanes", "any")
    lines += [f"{lib['fastq_id']},{lib['fastqs']},{lanes},{lib['feature_types']}" for lib in libraries]
    return lines


def main(smk):
    """Snakemake entry point."""
    redirect_to_log(smk.log[0])
    modalities = smk.params.modalities
    libs = smk.params.libraries
    libraries = []
    for m in modalities:
        sheet = load_samplesheet(smk.input[m])
        if libs[m] not in sheet:
            raise ValueError(
                f"{smk.wildcards.sample}: {m} library '{libs[m]}' has no FASTQs in {smk.input[m]} "
                f"(libraries found: {', '.join(sorted(sheet)) or 'none'})"
            )
        entry = sheet[libs[m]]
        libraries.append({**entry, "feature_types": FEATURE_TYPES[m]})
    with open(smk.output.csv, "w") as out:
        out.write("\n".join(config_lines(smk.config["cellranger"], libraries, modalities)) + "\n")
    log_message(f"Wrote {smk.output.csv} ({len(libraries)} libraries)", level="DONE")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
