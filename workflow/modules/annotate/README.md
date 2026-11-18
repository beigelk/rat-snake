# annotate

Cell metadata, then population by population: clusters, subsets, markers.

## Cell metadata (`metadata` rule)

Each dataset's `metadata` block (a spreadsheet, `.xlsx`/`.xls`/`.csv`/`.tsv`)
is read: `id_col` is matched to the sample names (Jaro-Winkler distance up to
`annotate.metadata.match_max_dist`; the matches are in
`results/annotate/metadata/sample_matches.tsv`), and `cols` are added to the
cells with `levels` as factor orders. With TCR data, scRepertoire
`combineExpression` adds the clonotype columns (`CTstrict`, `cloneSize`, ...).
Output: `results/annotate/metadata/cell_metadata.qs2` and `.tsv`.

## Populations (`annotate.populations`)

A population is a set of cells clustered at several resolutions. The first
is usually the whole integrated object; any other is cut from a parent's
clusters at one resolution and re-analysed (WNN and UMAP again):

```yaml
populations:
  - name: main
    resolutions:
      fine: {from: 0.05, to: 0.07, by: 0.01}     # or {values: [0.05, 0.1]}
    findallmarkers:
      assays: [RNA, ADT]
      threshold_sets:
        RNA: {default: {lfc_threshold: 0.25, padj_threshold: 0.05, pct1_threshold: 0.1}}
    plots:
      marker_sets:
        RNA:
          top5_default: {type: topn, threshold_set: default, top_n: 5}
          t_cells: {type: genelist, file: input/markers/t_cells.txt}
  - name: lymphocytes
    subset: {base_name: main, base_resolution: 0.07, base_clusters: 4}   # '0, 1, 8', [0, 1] or all
    resolutions:
      fine: {from: 0, to: 2, by: 0.2}
    wnn: {k_nn: 10}            # overrides annotate.subset_wnn for this population
```

Populations can subset subsets (main -> lymphocytes -> CD8s). A base
resolution the parent is not clustered at stops the run at start-up.

Per population, under `results/annotate/populations/{population}/`:

| Folder | Content |
|---|---|
| `cluster/` | `cluster__{population}_res.{res}.tsv`: barcode -> cluster |
| `subset/` | `subset__{population}.qs2` (subset populations) |
| `compile/` | `{population}_seurat.qs2` with every resolution's clusters and the cell metadata; `{population}_cluster_sizes.tsv` |
| `findallmarkers/` | FindAllMarkers per resolution and assay (`..._markers_raw.tsv`) |
| `filtermarkers/` | markers passing each threshold set (`..._{set}__markers_filtered.tsv`/`.xlsx`) |
| `clustree/` | clustree plot; `umaps/` has a cluster UMAP per resolution |

Cluster columns in the objects are `cluster__{population}__{graph}_res.{res}`.
Paths are built in `workflow/lib/rat_snake/populations.py`.

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `metadata` | rule | Build the cell metadata table: sample spreadsheets matched to samples, plus scRepertoire clonotypes. |
| `cluster` | rule | Cluster one population at one resolution (FindClusters on its WNN graph) and save the clusters. |
| `subset` | rule | Cut a population from its parent's clusters at one resolution and re-run WNN and UMAP. |
| `compile` | rule | One object per population with every resolution's clusters and the cell metadata; cells per cluster. |
| `findallmarkers` | rule | Positive markers of every cluster (FindAllMarkers) for one population, resolution and assay. |
| `filtermarkers` | rule | Keep the markers passing one threshold set (log2FC, adjusted p value, pct.1). |
| `clustree` | rule | clustree of a population's resolutions and a UMAP of its clusters at each one. |
<!-- RULES:END -->
