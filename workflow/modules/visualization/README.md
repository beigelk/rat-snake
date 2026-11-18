# visualization

- `metadata_plots`: UMAPs and QC violin/ridge plots of the integrated object
  coloured by sample metadata, from `annotate.metadata.plots` (`umap` and
  `feature_stat`, each with `defaults` and a list of `plots`). Figures go to
  `results/visualization/metadata/figures/`, listed in `metadata_plots.tsv`.
- `dotplot` and `markerplots`: for every population, resolution, assay and
  marker set (`plots.marker_sets` of each population), in
  `results/visualization/populations/{population}/`:
  - `dotplot/..__DotPlot.png` (Seurat) and `..__DotPlotClustered_k{n}.png`
    (scCustomize)
  - `markerplots/..__ComboPlot.png` (stacked violins and cluster UMAP above
    feature UMAPs), `..__StackViolin.png` and `..__FeatureDim.png`

A `topn` marker set takes the top `top_n` markers of each cluster from a
filtered marker table; a `genelist` set reads one feature per line. Features
not in the assay are logged and dropped; with none left, a placeholder figure
is written so the run goes on.

`visualization.reduction` (default `wnn.umap_harmony`) is the UMAP every plot
uses.

## Software

`env_r.yaml` (R 4.5, scCustomize, ComplexHeatmap). scplotter and plotthis
are only on CRAN, so `env_r.post.R` installs them after the env is built
(`envs/create_module_envs.sh` and the image build run it).

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `metadata_plots` | rule | UMAPs and QC violin/ridge plots of the integrated object (annotate.metadata.plots). |
| `dotplot` | rule | Dot plot and clustered dot plot of one marker set at one resolution. |
| `markerplots` | rule | Clusters UMAP, stacked violins and feature UMAPs of one marker set (one combined figure). |
<!-- RULES:END -->
