# preprocess

Per sample: ambient RNA, doublets, ADT normalisation, TCR clonotypes, then
one filtered RNA + ADT Seurat object. Every step writes a one-row summary,
joined into `results/preprocess/preprocess_summary.tsv` (cells kept at each
step, SoupX rho, doublets, dsb background drops, clonotypes).

Outputs are in `results/preprocess/{step}/{dataset}/` (`objects/`,
`tables/`, `figures/`).

| Step | What it does | Options |
|---|---|---|
| `soupx` | SCTransform, PCA and clusters on Cell Ranger's filtered cells; SoupX `autoEstCont` and `adjustCounts` | `preprocess.soupx` (`npcs`, `resolution`) |
| `doubletfinder` | Pre-filter cells, pK sweep, DoubletFinder with the 10x multiplet rate | `preprocess.doubletfinder` |
| `dsb` | Empty droplets from the raw matrix; dsb with isotype controls | `preprocess.dsb` |
| `tcr_clones` | scRepertoire `combineTCR` on Cell Ranger's filtered contigs | `preprocess.tcr` |
| `construct` | RNA from SoupX, doublets removed, RNA filter, ADT from dsb, IsoCount filter, ADT normalisation | `preprocess.construct`, `preprocess.mt_pattern` |

`preprocess.construct.adt.normalization` is `CLR` by default, as in the
original code (CLR on top of the dsb values); `dsb` keeps the dsb values.

TCR clonotypes run when `vdj` is in `modalities` and
`preprocess.tcr.enabled` is true; the clonotypes are added to the cells in
the annotate module's `metadata` rule (`combineExpression`).

## Software

`env_r.yaml` (R 4.4: SoupX, dsb, scRepertoire 2.2) and `env_doublet.yaml`
(R 4.5: DoubletFinder 2.0.6). dsb and scRepertoire 2.2 need R 4.4 on
conda-forge, DoubletFinder needs R 4.5, so DoubletFinder has its own env.

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `soupx` | rule | Remove ambient RNA with SoupX and build the sample's first Seurat object (RNA + ADT). |
| `doubletfinder` | rule | Pre-filter cells and call doublets with DoubletFinder (either of two runs calling a doublet counts). |
| `dsb` | rule | Normalise and denoise ADT counts with dsb, using empty drops as background and isotype controls. |
| `tcr_clones` | rule | Build the sample's TCR clonotypes from Cell Ranger's filtered contigs with scRepertoire combineTCR. |
| `construct` | rule | Build the sample's filtered RNA + ADT Seurat object: RNA QC, doublet removal, dsb values, ADT QC. |
| `preprocess_summary` | rule | Join every sample's SoupX, DoubletFinder, dsb, construct and TCR summaries into one QC table. |
<!-- RULES:END -->
