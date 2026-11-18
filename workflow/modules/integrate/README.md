# integrate

Merges every sample's object and integrates them.

- `merge`: one object (cells renamed `{sample}_{barcode}`).
- `harmony`: SCT and RNA PCA, ADT PCA (on `integrate.harmony.adt_markers_file`
  features, or every non-isotype ADT), Harmony on both by
  `integrate.harmony.group_by_vars`, then WNN (`integrate.wnn`) and the WNN
  UMAP. Output: `results/integrate/harmony/harmonized_seurat.qs2`, the root
  object every population starts from.

`integrate.wnn.prune_snn` takes a number or a fraction such as `"1/15"`;
`rna_dims`/`adt_dims` default to every Harmony dimension.

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `merge` | rule | Merge every usable sample's Seurat object; cell names become {sample}_{barcode}. |
| `harmony` | rule | Harmony on the RNA and ADT PCAs, then WNN graphs (wsnn_harmony) and the wnn.umap_harmony UMAP. |
<!-- RULES:END -->
