"""Unit tests for workflow/lib/rat_snake (samples and populations)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "lib"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "core", "python"))

from rat_snake.populations import Populations, expand_resolutions, format_res, parse_clusters  # noqa: E402
from rat_snake.samples import (  # noqa: E402
    build_samples,
    cellranger_files,
    fastq_libraries,
    read_sample_table,
    resolve_run_dir,
    shared_tokens,
)
from snake_core.paths import ProjectPaths  # noqa: E402


# -- samples -------------------------------------------------------------------
def test_shared_tokens_drops_only_sequencer_numbers():
    assert shared_tokens(["25_13_1205_S25", "25_13_1205_S41"]) == "25_13_1205"
    # The old code dropped every part starting with S
    assert shared_tokens(["Spleen_7_S1", "Spleen_7_S2"]) == "Spleen_7"


def write_table(path, rows, header="gex\tcsp\tvdj"):
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return str(path)


def test_read_sample_table(tmp_path):
    t = write_table(tmp_path / "s.tsv", ["A_1_S1\tA_1_S2\tA_1_S3", "B_S4\tB_S5\t"])
    s = read_sample_table(t)
    assert list(s) == ["A_1", "B"]
    assert s["B"] == {"gex": "B_S4", "csp": "B_S5"}


def test_read_sample_table_sample_column_and_duplicates(tmp_path):
    t = write_table(tmp_path / "s.tsv", ["x_S1\tx_S2\tx_S3\tmine"], header="gex\tcsp\tvdj\tsample")
    assert list(read_sample_table(t)) == ["mine"]
    t = write_table(tmp_path / "d.tsv", ["A_S1\tA_S2\tA_S3", "A_S4\tA_S5\tA_S6"])
    with pytest.raises(ValueError, match="used twice"):
        read_sample_table(t)
    t = write_table(tmp_path / "m.tsv", ["A_S1\tA_S2"], header="gex\tcsp")
    with pytest.raises(ValueError, match="missing column"):
        read_sample_table(t)


def make_run(root, sample, vdj=True):
    outs = root / sample / "outs"
    (outs / "multi" / "count" / "raw_feature_bc_matrix").mkdir(parents=True)
    per = outs / "per_sample_outs" / f"other_{sample}"
    (per / "count" / "sample_filtered_feature_bc_matrix").mkdir(parents=True)
    (per / "metrics_summary.csv").write_text("x\n")
    if vdj:
        (per / "vdj_t").mkdir()
        (per / "vdj_t" / "filtered_contig_annotations.csv").write_text("x\n")
    return root / sample


def test_resolve_run_dir_and_files(tmp_path):
    run = make_run(tmp_path, "A_1")
    assert resolve_run_dir(run) == str(run)
    assert resolve_run_dir(run / "outs") == str(run)
    with pytest.raises(FileNotFoundError):
        resolve_run_dir(tmp_path / "nope")
    # The only per_sample_outs folder is used when none is named after the sample
    files = cellranger_files(str(run), "A_1")
    assert files["metrics"].endswith("other_A_1/metrics_summary.csv")


def write_fastqs(folder, libs):
    """Empty R1/R2 FASTQs named like the sequencer's: {TAG}_{library}_L001_R1_001.fastq.gz."""
    folder.mkdir(exist_ok=True)
    for tag, lib in libs:
        for r in ("R1", "R2"):
            (folder / f"{tag}_{lib}_L001_{r}_001.fastq.gz").write_text("")


def config_for(tmp_path, existing=None, modalities=("gex", "csp", "vdj")):
    t = write_table(tmp_path / "s.tsv", ["A_1_S1\tA_1_S2\tA_1_S3", "B_S4\tB_S5\t"])
    fq = tmp_path / "fastq"
    write_fastqs(fq, [("GEX", "A_1_S1"), ("CSP", "A_1_S2"), ("VDJ", "A_1_S3"), ("GEX", "B_S4"), ("CSP", "B_S5")])
    return {
        "modalities": list(modalities),
        "datasets": {"d1": {"sample_ids": t, "input_dirs": {"gex": str(fq), "csp": str(fq), "vdj": str(fq)}},
                     "off": {"sample_ids": "x", "enabled": False}},
        "cellranger": {"existing_outs": existing or {}},
    }


def test_build_samples_statuses(tmp_path):
    s = build_samples(config_for(tmp_path))
    assert list(s) == ["d1"]
    assert s["d1"]["A_1"]["status"] == "run"
    assert s["d1"]["B"]["status"] == "excluded"
    # Without vdj in modalities, B is usable
    s = build_samples(config_for(tmp_path, modalities=("gex", "csp")))
    assert s["d1"]["B"]["status"] == "run"


def test_samples_without_fastqs_are_excluded(tmp_path):
    cfg = config_for(tmp_path)
    # Only R1 for A_1's VDJ library: not enough
    (tmp_path / "fastq" / "VDJ_A_1_S3_L001_R2_001.fastq.gz").unlink()
    s = build_samples(cfg)
    assert s["d1"]["A_1"]["status"] == "excluded"
    assert "vdj (A_1_S3)" in s["d1"]["A_1"]["reason"]
    # A finished Cell Ranger run needs no FASTQs
    root = tmp_path / "runs"
    make_run(root, "A_1")
    cfg["cellranger"]["existing_outs"] = {"root": str(root)}
    assert build_samples(cfg)["d1"]["A_1"]["status"] == "existing"


def test_fastq_libraries_match_samplesheet(tmp_path):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "modules", "setup", "scripts"))
    from make_samplesheet import scan

    fq = tmp_path / "fq"
    write_fastqs(fq, [("GEX", "25_13_1205_S25"), ("GEX", "Spleen_7_S1"), ("CSP", "25_13_1205_S41")])
    (fq / "GEX-only_r1_S9_L002_R1_001.fastq.gz").write_text("")
    data, _, _ = scan([str(fq)], "fastq.gz", "GEX")
    paired = {lib for lib, lanes in data.items() if any({"R1", "R2"} <= set(r) for r in lanes.values())}
    assert fastq_libraries([str(fq)], "fastq.gz", "gex") == paired == {"25_13_1205_S25", "Spleen_7_S1"}


def test_build_samples_existing_outs(tmp_path):
    root = tmp_path / "runs"
    make_run(root / "d1", "A_1")
    s = build_samples(config_for(tmp_path, {"root": str(root)}))
    assert s["d1"]["A_1"]["status"] == "existing"
    assert s["d1"]["A_1"]["run_dir"] == str(root / "d1" / "A_1")
    # A run named explicitly but missing contigs fails at start-up
    bad = make_run(tmp_path / "bad", "A_1", vdj=False)
    with pytest.raises(ValueError, match="missing"):
        build_samples(config_for(tmp_path, {"samples": {"d1/A_1": str(bad)}}))


def test_missing_sample_table(tmp_path):
    cfg = config_for(tmp_path)
    cfg["datasets"]["d1"]["sample_ids"] = str(tmp_path / "none.tsv")
    with pytest.raises(FileNotFoundError, match="datasets.d1.sample_ids"):
        build_samples(cfg)


# -- populations -----------------------------------------------------------------
def test_format_res_matches_r():
    assert format_res(0.0) == "0"
    assert format_res(1.0) == "1"
    assert format_res(0.10) == "0.1"
    assert format_res(0.07) == "0.07"
    assert format_res("{res}") == "{res}"


def test_expand_resolutions_and_clusters():
    assert expand_resolutions({"from": 0, "to": 1, "by": 0.2}) == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    assert expand_resolutions({"values": [0.07, 0.1]}) == [0.07, 0.1]
    with pytest.raises(KeyError):
        expand_resolutions({"from": 0})
    assert parse_clusters("0, 1, 8") == [0, 1, 8]
    assert parse_clusters(4) == [4]
    assert parse_clusters([0, "1"]) == [0, 1]
    assert parse_clusters("All") == "all"


def pops_config(base_res=0.07):
    fam = {"assays": ["RNA"], "threshold_sets": {"RNA": {"default": {}, "strict": {}}}}
    plots = {"marker_sets": {"RNA": {"top5": {"type": "topn", "threshold_set": "strict", "top_n": 5},
                                     "list": {"type": "genelist", "file": "m.txt"}}}}
    return {"annotate": {"populations": [
        {"name": "main", "resolutions": {"a": {"from": 0.05, "to": 0.07, "by": 0.01}, "b": {"values": [1]}},
         "findallmarkers": fam, "plots": plots},
        {"name": "sub", "graph_name": "g2", "subset": {"base_name": "main", "base_resolution": base_res,
                                                        "base_clusters": "0, 2"},
         "resolutions": {"a": {"values": [0.5]}}, "findallmarkers": fam},
    ]}}


def test_populations_plan_and_paths(tmp_path):
    P = ProjectPaths({"project": {"name": "p", "version": "v1", "output_root": str(tmp_path)}})
    pops = Populations(pops_config(), P, "root.qs2")
    assert pops.names == ["main", "sub"]
    assert pops.resolutions("main") == [0.05, 0.06, 0.07, 1.0]
    assert pops.is_root("main") and not pops.is_root("sub")
    assert pops.input_object("main") == "root.qs2"
    assert pops.input_object("sub").endswith("annotate/populations/sub/subset/subset__sub.qs2")
    assert pops.subset_spec("sub") == {"base_name": "main", "base_graph_name": "wsnn_harmony",
                                       "base_resolution": "0.07", "base_clusters": [0, 2]}
    assert pops.cluster_column("sub", 0.5) == "cluster__sub__g2_res.0.5"
    assert pops.cluster_tsv("main", 1.0).endswith("cluster__main_res.1.tsv")
    assert len(pops.marker_jobs()) == 5
    assert len(pops.filter_jobs()) == 10
    assert len(pops.plot_jobs()) == 8
    assert pops.plot_features_input("main", 1.0, "RNA", "top5").endswith("p_v1__main__RNA_res.1__strict__markers_filtered.tsv")
    assert pops.plot_features_input("main", 1.0, "RNA", "list") == os.path.abspath("m.txt")
    assert pops.wnn("sub", {"k_nn": 10, "rna_dims": 15}) == {"k_nn": 10, "rna_dims": 15}


def test_populations_bad_base_resolution(tmp_path):
    P = ProjectPaths({"project": {"name": "p", "version": "v1", "output_root": str(tmp_path)}})
    with pytest.raises(ValueError, match="not clustered at"):
        Populations(pops_config(base_res=0.3), P, "root.qs2")
