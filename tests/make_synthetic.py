"""Make a small synthetic CITE-seq + TCR project for testing rat-snake.

Writes, under the output folder (default tests/data):

    fastq/{dataset}/          empty FASTQs named like the sequencer's, so the
                              setup and cellranger modules can be dry-run
    cellranger/{sample}/outs/ finished Cell Ranger multi runs (raw and
                              filtered matrices, contigs, metrics), used
                              through cellranger.existing_outs
    sample_ids_{dataset}.tsv  sample tables (one sample lacks a VDJ library,
                              so it is excluded)
    metadata_{dataset}.csv    sample metadata sheets
    markers/*.txt             gene lists for genelist marker sets
    config_existing.yaml      config that reads the finished runs
    config_run.yaml           config that would run Cell Ranger (dry-run only)

The cells come from four cell types with their own genes and proteins,
plus empty droplets with low RNA and some ADT counts, so SoupX, DoubletFinder,
dsb, clustering and markers all have something to find.

Usage:
    python tests/make_synthetic.py [--out tests/data] [--cells 400] [--seed 1]
"""

import argparse
import csv
import gzip
import os

import numpy as np
import yaml

N_GENES = 600
N_MT = 13
CELL_TYPES = 4
GENES_PER_TYPE = 10
ADTS = [f"CD{i}" for i in (3, 4, 8, 14, 19, 56, 45, 11, 16, 25, 27, 28, 38, 44, 69, 127)] + [
    "Mouse-IgG1-Isotype",
    "Mouse-IgG2a-Isotype",
    "Rat-IgG2b-Isotype",
]

# dataset -> [(sample, gex, csp, vdj)]; S numbers are the sequencer's
SAMPLES = {
    "d1": [
        ("25_13_1205", "25_13_1205_S25", "25_13_1205_S41", "25_13_1205_S57"),
        ("25_13_1206", "25_13_1206_S26", "25_13_1206_S42", "25_13_1206_S58"),
        ("25_13_1207", "25_13_1207_S27", "25_13_1207_S43", ""),  # no VDJ: excluded
    ],
    "d2": [
        ("Spleen_7", "Spleen_7_S1", "Spleen_7_S2", "Spleen_7_S3"),
    ],
}


def gene_names():
    """Gene symbols: MT- genes, cell type genes and background genes."""
    genes = [f"MT-G{i}" for i in range(N_MT)]
    for t in range(CELL_TYPES):
        genes += [f"TYPE{t}G{i}" for i in range(GENES_PER_TYPE)]
    genes += [f"GENE{i}" for i in range(N_GENES - len(genes))]
    return genes


def simulate(rng, n_cells, n_empty):
    """Counts of cells and empty droplets.

    Returns:
        (rna, adt, cell_types): genes x droplets and ADTs x droplets count
        matrices (cells first) and each cell's type.
    """
    genes = gene_names()
    # Expression profiles are shared by every sample (fixed seed), so samples
    # differ only by noise and integrate cleanly
    shared = np.random.default_rng(0)
    base = shared.gamma(2.0, 1.0, len(genes))
    base[:N_MT] *= 3
    base[N_MT : N_MT + CELL_TYPES * GENES_PER_TYPE] = 0.01  # cell type genes: off outside their type
    types = rng.integers(0, CELL_TYPES, n_cells)
    rna = np.zeros((len(genes), n_cells + n_empty), dtype=np.int64)
    adt = np.zeros((len(ADTS), n_cells + n_empty), dtype=np.int64)
    adt_base = shared.gamma(2.0, 1.0, len(ADTS))
    adt_base[-3:] = 0.5  # isotype controls stay low
    for j in range(n_cells):
        prof = base.copy()
        lo = N_MT + types[j] * GENES_PER_TYPE
        prof[lo : lo + GENES_PER_TYPE] = 30.0
        prof /= prof.sum()
        rna[:, j] = rng.poisson(prof * rng.integers(3000, 6000))
        aprof = adt_base * rng.lognormal(0, 0.5, len(ADTS))
        aprof[types[j] * 4 : types[j] * 4 + 4] *= 30
        aprof /= aprof.sum()
        adt[:, j] = rng.poisson(aprof * rng.integers(1500, 3000))
    # Ambient RNA is the pooled cells' RNA, so cell type genes are in the soup
    soup = rna[:, :n_cells].sum(axis=1) / rna[:, :n_cells].sum()
    for j in range(n_cells, n_cells + n_empty):
        rna[:, j] = rng.poisson(soup * rng.integers(5, 150))
        adt[:, j] = rng.poisson(adt_base / adt_base.sum() * rng.integers(40, 600))
    # A little soup in every cell, for SoupX to estimate
    for j in range(n_cells):
        rna[:, j] += rng.poisson(soup * 15)
    return rna, adt, types


def write_matrix(folder, rna, adt, barcodes):
    """A Cell Ranger feature-barcode matrix folder (gzipped MTX)."""
    os.makedirs(folder, exist_ok=True)
    genes = gene_names()
    with gzip.open(os.path.join(folder, "features.tsv.gz"), "wt") as fh:
        for i, g in enumerate(genes):
            fh.write(f"ENSG{i:08d}\t{g}\tGene Expression\n")
        for a in ADTS:
            fh.write(f"{a}\t{a}\tAntibody Capture\n")
    with gzip.open(os.path.join(folder, "barcodes.tsv.gz"), "wt") as fh:
        fh.writelines(f"{b}\n" for b in barcodes)
    mat = np.vstack([rna, adt])
    rows, cols = np.nonzero(mat)
    with gzip.open(os.path.join(folder, "matrix.mtx.gz"), "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        fh.write(f"{mat.shape[0]} {mat.shape[1]} {len(rows)}\n")
        for r, c in zip(rows, cols):
            fh.write(f"{r + 1} {c + 1} {mat[r, c]}\n")


def write_contigs(path, rng, barcodes, types):
    """filtered_contig_annotations.csv: TRA and TRB chains for most T cells (types 0 and 1)."""
    cols = ["barcode", "is_cell", "contig_id", "high_confidence", "length", "chain", "v_gene", "d_gene",
            "j_gene", "c_gene", "full_length", "productive", "cdr3", "cdr3_nt", "reads", "umis",
            "raw_clonotype_id", "raw_consensus_id"]
    aa = "ACDEFGHIKLMNPQRSTVWY"
    clones = []
    for _ in range(40):
        cdr3 = ("CA" + "".join(rng.choice(list(aa), 9)) + "F", "CAS" + "".join(rng.choice(list(aa), 9)) + "F")
        nts = tuple("".join(rng.choice(list("ACGT"), len(c) * 3)) for c in cdr3)
        clones.append((cdr3, nts))
    weights = np.r_[np.full(5, 8.0), np.ones(35)]
    weights /= weights.sum()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for b, t in zip(barcodes, types):
            if t > 1 or rng.random() < 0.15:
                continue
            k = rng.choice(len(clones), p=weights)
            (cdr3_a, cdr3_b), (nt_a, nt_b) = clones[k]
            for chain, cdr3, nt, v, j, c in (("TRA", cdr3_a, nt_a, f"TRAV{k % 20 + 1}", "TRAJ1", "TRAC"),
                                             ("TRB", cdr3_b, nt_b, f"TRBV{k % 25 + 1}", "TRBJ1-1", "TRBC1")):
                w.writerow([b, "true", f"{b}_contig_{chain}", "true", 500, chain, v, "", j, c, "true", "true",
                            cdr3, nt, 1000, 10, f"clonotype{k + 1}", f"clonotype{k + 1}_consensus_1"])


def write_metrics(path, n_cells):
    """metrics_summary.csv in Cell Ranger multi's long format."""
    rows = [
        ("Cells", "Gene Expression", "", "", "Cells", str(n_cells)),
        ("Cells", "Gene Expression", "", "", "Median reads per cell", "20,000"),
        ("Cells", "Gene Expression", "", "", "Median genes per cell", "900"),
        ("Cells", "Antibody Capture", "", "", "Median UMI counts per cell", "2,000"),
        ("Cells", "VDJ T", "", "", "Cells with productive V-J spanning (TRA, TRB) pair", "40%"),
        ("Library", "Gene Expression", "Physical library ID", "GEX_1", "Sequencing saturation", "60%"),
    ]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Category", "Library Type", "Grouped By", "Group Name", "Metric Name", "Metric Value"])
        w.writerows(rows)


def make_run(run_dir, sample, rng, n_cells):
    """One finished Cell Ranger multi run folder."""
    n_empty = n_cells * 2
    rna, adt, types = simulate(rng, n_cells, n_empty)
    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n_cells + n_empty)]
    outs = os.path.join(run_dir, "outs")
    write_matrix(os.path.join(outs, "multi", "count", "raw_feature_bc_matrix"), rna, adt, barcodes)
    per = os.path.join(outs, "per_sample_outs", sample)
    write_matrix(os.path.join(per, "count", "sample_filtered_feature_bc_matrix"),
                 rna[:, :n_cells], adt[:, :n_cells], barcodes[:n_cells])
    os.makedirs(os.path.join(per, "vdj_t"), exist_ok=True)
    write_contigs(os.path.join(per, "vdj_t", "filtered_contig_annotations.csv"), rng, barcodes[:n_cells], types)
    write_metrics(os.path.join(per, "metrics_summary.csv"), n_cells)
    with open(os.path.join(per, "web_summary.html"), "w") as fh:
        fh.write(f"<html><body>{sample}</body></html>\n")


def populations():
    """A root population and one subset, with topn and genelist marker sets."""
    fam = {"assays": ["RNA", "ADT"], "threshold_sets": {
        "RNA": {"default": {"lfc_threshold": 0.25, "padj_threshold": 0.05, "pct1_threshold": 0.1}},
        "ADT": {"default": {"lfc_threshold": 0.5, "padj_threshold": 0.05, "pct1_threshold": 0.5}}}}
    sets = {"RNA": {"top5_default": {"type": "topn", "threshold_set": "default", "top_n": 5},
                    "type0": {"type": "genelist", "file": "tests/data/markers/type0.txt"}},
            "ADT": {"top5_default": {"type": "topn", "threshold_set": "default", "top_n": 5}}}
    return [
        {"name": "main", "resolutions": {"fine": {"values": [0.1, 0.3]}}, "findallmarkers": fam,
         "plots": {"marker_sets": sets}},
        {"name": "tcells", "subset": {"base_name": "main", "base_resolution": 0.1, "base_clusters": "0, 1"},
         "resolutions": {"fine": {"from": 0.2, "to": 0.4, "by": 0.2}}, "findallmarkers": fam,
         "plots": {"marker_sets": sets}, "wnn": {"rna_dims": 10, "adt_dims": 8}},
    ]


def config(out, existing):
    """Test config: every dataset, the synthetic references and populations."""
    cfg = {
        "project": {"name": "synthetic", "version": "v0", "output_root": os.path.join(out, "output")},
        "software": {"mode": "none"},
        "r_options": {"future_globals_gb": 4},
        "datasets": {},
        "cellranger": {"gex_reference": os.path.join(out, "refs", "gex"),
                       "csp_reference": os.path.join(out, "refs", "feature_ref.csv"),
                       "vdj_reference": os.path.join(out, "refs", "vdj"),
                       "path": "cellranger", "threads": 2, "mem_gb": 4},
        "preprocess": {
                        "doubletfinder": {"min_feat_cutoff": 100, "workers": 1},
            "construct": {"rna": {"rna_min_feat_cutoff": 100}, "adt": {"adt_feat_min": 5}, "npcs": 10},
        },
        "integrate": {"harmony": {"rna_npcs": 10, "adt_npcs": 8}, "wnn": {"k_nn": 15}},
        "annotate": {
            "subset_wnn": {"k_nn": 10, "prune_snn": "1/30", "rna_dims": 10, "adt_dims": 8},
            "metadata": {"plots": {"umap": {"defaults": {"width": 8, "height": 6},
                                            "plots": [{"group_var": "Dataset"}, {"group_var": "Condition"}]},
                                   "feature_stat": {"features": ["nCount_RNA", "percent_mt"],
                                                    "defaults": {"width": 8, "height": 5},
                                                    "plots": [{"plot_type": "violin", "group_var": "Sample_ID"}]}}},
            "populations": populations(),
        },
        "report": {"resolutions": {"main": [0.1]}},
    }
    for d, rows in SAMPLES.items():
        cfg["datasets"][d] = {
            "sample_ids": os.path.join(out, f"sample_ids_{d}.tsv"),
            "input_dirs": {m: os.path.join(out, "fastq", d) for m in ("gex", "csp", "vdj")},
            "metadata": {"file": os.path.join(out, f"metadata_{d}.csv"), "id_col": "Sample name",
                         "cols": ["Sample name", "Condition"]},
        }
    if existing:
        cfg["cellranger"]["existing_outs"] = {"root": os.path.join(out, "cellranger")}
    return cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="tests/data")
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    out = args.out
    rng = np.random.default_rng(args.seed)

    for d, rows in SAMPLES.items():
        fq = os.path.join(out, "fastq", d)
        os.makedirs(fq, exist_ok=True)
        with open(os.path.join(out, f"sample_ids_{d}.tsv"), "w") as fh:
            fh.write("gex\tcsp\tvdj\n")
            for sample, gex, csp, vdj in rows:
                fh.write(f"{gex}\t{csp}\t{vdj}\n")
                for tag, lib in (("GEX", gex), ("CSP", csp), ("VDJ", vdj)):
                    if lib:
                        for r in ("R1", "R2"):
                            open(os.path.join(fq, f"{tag}_{lib}_L001_{r}_001.fastq.gz"), "w").close()
                if vdj:
                    make_run(os.path.join(out, "cellranger", sample), sample, rng, args.cells)
        with open(os.path.join(out, f"metadata_{d}.csv"), "w") as fh:
            fh.write("Sample name,Condition\n")
            for i, (sample, *_rest) in enumerate(rows):
                fh.write(f"{sample},{'Treated' if i % 2 else 'Control'}\n")

    os.makedirs(os.path.join(out, "markers"), exist_ok=True)
    with open(os.path.join(out, "markers", "type0.txt"), "w") as fh:
        fh.write("\n".join(f"TYPE0G{i}" for i in range(6)) + "\nNOT_A_GENE\n")
    for name in ("gex", "vdj"):
        os.makedirs(os.path.join(out, "refs", name), exist_ok=True)
    with open(os.path.join(out, "refs", "feature_ref.csv"), "w") as fh:
        fh.write("id,name,read,pattern,sequence,feature_type\n")
        fh.writelines(f"{a},{a},R2,5PNNNNNNNNNN(BC),ACGTACGTACGT,Antibody Capture\n" for a in ADTS)

    for name, existing in (("config_existing.yaml", True), ("config_run.yaml", False)):
        with open(os.path.join(out, name), "w") as fh:
            yaml.safe_dump(config(out, existing), fh, sort_keys=False)
    print(f"Synthetic project written to {out}")


if __name__ == "__main__":
    main()
