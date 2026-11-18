# cellranger

Runs Cell Ranger multi for every sample with status `run`, or uses finished
runs, then collects their metrics and web summaries.

- `results/cellranger/configs/{dataset}/{sample}_config.csv`: the multi
  config (references, the sample's GEX, CSP and VDJ libraries).
- `results/cellranger/multi/{dataset}/{sample}/`: the run (`outs/`).
- `results/cellranger/metrics/cellranger_metrics.tsv`: one row per sample,
  one column per metric (`_long.tsv`: sample, metric, value).
- `results/cellranger/web_summaries/{project}_CellRanger_web_summaries.zip`.

## Running Cell Ranger

Cell Ranger is not a conda package. Name it in the config, as a path or as an
environment module:

```yaml
cellranger:
  path: /opt/cellranger-9.0.1/bin/cellranger   # or leave null and:
software:
  envmodules:
    cellranger: [CellRanger/9.0.1]             # adds --use-envmodules
```

`cellranger.threads` and `cellranger.mem_gb` set `--localcores`/`--localmem`
and the job's resources. With `cellranger.scratch_dir` (e.g.
`/scr1/users/$USER`), each job runs in its own folder there and copies the
results back; `SC_MULTI_CS` is deleted unless `keep_sc_multi_cs: true`.

## Using finished runs (`cellranger.existing_outs`)

To skip Cell Ranger, for debugging or to reuse earlier runs, point the
pipeline at finished runs. A sample found here is never re-run; its outputs
are read where they are (nothing is copied).

```yaml
cellranger:
  existing_outs:
    root: /data/old_project/cellranger      # searched for {root}/{dataset}/{sample} and {root}/{sample}
    samples:                                # or name runs one by one (overrides root)
      25_13_1205: /data/runs/25_13_1205             # run folder (holding outs/) ...
      dataset2/Spleen_7: /data/runs/Spleen_7/outs   # ... or its outs/; key by sample or dataset/sample
```

At start-up each run is checked for the raw matrix, the filtered matrix,
`metrics_summary.csv` and (with TCR) `filtered_contig_annotations.csv`; a
missing file stops the run with its path. `per_sample_outs/{sample}` is used,
or the only folder under `per_sample_outs` when the run used another name.
`results/setup/sample_status.tsv` lists which run each sample uses.

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `make_cellranger_config` | rule | Write one sample's Cell Ranger multi config from the samplesheets. |
| `cellranger_multi` | rule | Run Cell Ranger multi for one sample. |
| `cellranger_metrics` | rule | Combine every sample's Cell Ranger metrics_summary.csv into one table (sample x metric). |
| `collect_web_summaries` | rule | Zip every sample's Cell Ranger web_summary.html, renamed to {dataset}__{sample}_web_summary.html. |
<!-- RULES:END -->
