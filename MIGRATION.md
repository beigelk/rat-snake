# Moving from the dev layout

What changed when rat-snake was restructured onto snake-core (branch
`restructure`, from `dev`). The analysis is the same; file locations, config
keys and some output names changed, and wiring bugs were fixed.

## Running

| Before | Now |
|---|---|
| `bash rat-snake.sh -c config.yaml -r run` | same, plus `-e local|slurm`, `-s` (use the active env) and `-- <snakemake args>` |
| `snakemake/envs/conda/conda.sh` | `envs/conda_setup.sh` (Snakemake env) and `envs/create_module_envs.sh` (module envs) |
| one R module (`R/4.5.3`) and `resources.r_library` for every rule | per-module conda envs, images or `none` (`software.mode`); `software.r_library` still adds a library |
| Cell Ranger only through SLURM | runs locally or through the slurm profile; `cellranger.path` or `software.envmodules.cellranger` |

## Files

| Before | Now |
|---|---|
| `snakemake/Snakefile` | `workflow/Snakefile` |
| `snakemake/config/project/project_config_template.yaml` | `config/config.yaml` (+ `config/config.schema.yaml`) |
| `snakemake/profiles/` | `profiles/` |
| `snakemake/workflow/rules/<step>/*.smk` | `workflow/modules/<module>/<module>.smk` |
| `snakemake/workflow/scripts/<step>/*.R` | `workflow/modules/<module>/scripts/*.R` |
| `scripts/utils/{file_io,plot,utils}.R` | `core/R/{io,logging,snake_init}.R` and `workflow/lib/R/` |
| `preprocess/tcr/clones.R` | `preprocess/scripts/tcr_clones.R` |
| `visualization/combo_plot.R` | `visualization/scripts/markerplots.R` |
| `visualization/metadata.R` | plotting: `visualization/scripts/metadata_plots.R`; metadata table: `annotate/scripts/metadata.R` |
| `visualization/clustree.R` | `annotate/scripts/clustree.R` |
| `utils/cellranger_metrics_summary.R` | `cellranger/scripts/cellranger_metrics.R` |
| `snakemake/config/base.yaml`, `envs/docker/Dockerfile` | removed (`core/docker/conda.Dockerfile` builds the images) |

## Config keys

| Before | Now |
|---|---|
| `project_info: {name, version}` | `project: {name, version, output_root}` |
| `resources.r_library` | `software.r_library` |
| `datasets.<d>.metadata.input: {...}` | `datasets.<d>.metadata: {...}` (same fields; CSV/TSV also read) |
| `preprocess.rna.soupx` | `preprocess.soupx` |
| `preprocess.rna.doubletfinder` | `preprocess.doubletfinder` (+ `pN`, `workers`) |
| `preprocess.adt.dsb` | `preprocess.dsb` (thresholds that were in the script) |
| `preprocess.tcr.clones` | `preprocess.tcr` (scRepertoire options) |
| `preprocess.multimodal.construct` | `preprocess.construct` (+ `adt.normalization`, `npcs`, `cluster_resolution`) |
| hard-coded mitochondrial pattern | `preprocess.mt_pattern` (`^MT-`; mouse `^mt-`) |
| hard-coded WNN settings in harmony.R and subset.R | `integrate.wnn`, `annotate.subset_wnn`, a population's `wnn` |
| — | `seed`, `debug`, `modules`, `modalities`, `r_options`, `cellranger.existing_outs`, `cellranger.scratch_dir`, `visualization`, `report` |

`annotate.populations` and `annotate.metadata.plots` keep their format.

## Outputs

Everything is under `output/{name}_{version}/` (was `projects/...`):
`results/{module}/...`, `logs/`, `benchmarks/`, `reports/`, `debug/`.

| Before | Now |
|---|---|
| per-sample folders under each step | `results/preprocess/{step}/{dataset}/{objects,tables,figures}/` |
| `.done` files for plots | the plots themselves are the rule outputs |
| `sample.nameID` metadata column | `Sample_ID` (plus `Dataset`) |
| cluster files with the graph name | `results/annotate/populations/{population}/...` (file names leave out the graph name; cluster columns keep it: `cluster__{population}__{graph}_res.{res}`) |
| — | `results/setup/sample_status.tsv`, `results/preprocess/preprocess_summary.tsv`, `results/annotate/metadata/sample_matches.tsv`, `reports/report.{html,docx}` |

## Fixed

- Cell Ranger outputs were not wired to the preprocessing rules, and `rule
  all` was commented out, so the pipeline did not run end to end.
- The samplesheet checkpoint ran at parse time; sample status is now read
  from the sample tables at start-up, and samples missing a library are
  listed instead of failing.
- Sample names dropped every name part starting with "S" (e.g. "Spleen");
  now only sequencer numbers (`S` + digits) are dropped.
- The TCR rule and the marker plot rule referred to missing scripts and
  wildcards; dotplot.R used undefined `assay_name`, `donefile` and
  `r_library`; combo_plot.R wrote to a fixed `V5` folder and used
  `row2_t_margin` before setting it.
- Sample metadata joins dropped the cell names (all metadata NA); compile
  never added the metadata to the population objects; metadata.R saved its
  table as qs2, then overwrote it with RDS.
- `percent_mt` was truncated to an integer; `all.equal()` inside `if()`;
  `Sample_ID` was used but never created.
- The clustered dot plot failed when a marker set had fewer features than
  clusters; the significant-PC count was NA (and the next step failed) when
  either cut-off rule found nothing.
- Cell Ranger config paths are absolute (Cell Ranger runs in its own folder).

## Kept on purpose

- ADT normalisation stays CLR on top of the dsb values
  (`preprocess.construct.adt.normalization: CLR`); `dsb` keeps the dsb values.
- Filtered marker tables keep FindAllMarkers' order, so `topn` marker sets
  pick the same genes as before.
- scRepertoire defaults (`proportion: true`, strict clone call).
