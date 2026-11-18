"""Checks on a finished run of the synthetic project (tests/run_synthetic.sh).

The run folder defaults to tests/data/output/synthetic_v0; set RAT_OUTPUT to
check another one.
"""

import os

import pandas as pd
import pytest

OUT = os.environ.get("RAT_OUTPUT", os.path.join(os.path.dirname(__file__), "data", "output", "synthetic_v0"))
RES = os.path.join(OUT, "results")

pytestmark = pytest.mark.skipif(not os.path.isdir(RES), reason=f"no synthetic run at {OUT}")


def tsv(*parts):
    return pd.read_csv(os.path.join(RES, *parts), sep="\t")


def test_sample_status():
    st = tsv("setup", "sample_status.tsv").set_index("sample")
    assert st.loc["25_13_1207", "status"] == "excluded"
    assert set(st.drop("25_13_1207")["status"]) == {"existing"}


def test_preprocessing_keeps_cells_and_finds_clones():
    qc = tsv("preprocess", "preprocess_summary.tsv")
    assert list(qc["sample"]) == ["25_13_1205", "25_13_1206", "Spleen_7"]
    assert (qc["cells_final"] > 300).all()
    assert (qc["background_drops"] > 100).all()
    # Clones are shared by several cells, so there are fewer clonotypes than T cells
    assert (qc["clonotypes"] < qc["tcr_cells"]).all()


def test_four_cell_types_found():
    sizes = tsv("annotate", "populations", "main", "compile", "main_cluster_sizes.tsv")
    assert sizes.groupby("column")["cluster"].nunique().max() == 4


def test_markers_are_the_cell_type_genes():
    f = os.path.join(RES, "annotate", "populations", "main", "filtermarkers",
                     "synthetic_v0__main__RNA_res.0.1__default__markers_filtered.tsv")
    markers = pd.read_csv(f, sep="\t")
    top = markers.groupby("cluster").head(5)
    assert top["gene"].str.startswith("TYPE").all()
    # Each cluster's top markers come from one cell type
    assert (top.groupby("cluster")["gene"].apply(lambda g: g.str[:5].nunique()) == 1).all()


def test_subset_population():
    sizes = tsv("annotate", "populations", "tcells", "compile", "tcells_cluster_sizes.tsv")
    assert sizes.groupby("column")["cells"].sum().iloc[0] > 400


def test_metadata_matched():
    m = tsv("annotate", "metadata", "sample_matches.tsv")
    assert len(m) == 3 and (m["distance"] == 0).all()


def test_plots_and_report():
    assert os.path.exists(os.path.join(RES, "visualization", "populations", "main", "dotplot",
                                       "synthetic_v0__main__RNA_res.0.1__type0__DotPlot.png"))
    for f in ("report.html", "report.docx", "report_notion.docx"):
        assert os.path.getsize(os.path.join(OUT, "reports", f)) > 10_000
