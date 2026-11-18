# setup

Records what a run used and lists the FASTQs Cell Ranger will read.

- `config/config.yaml` in the output folder is the config after defaults are filled in.
- `results/setup/sample_status.tsv` lists every sample in every dataset's
  `sample_ids` table: its libraries, whether Cell Ranger runs for it (`run`),
  a finished run is used (`existing`, from `cellranger.existing_outs`) or it
  is left out (`excluded`, with the reason).
- `results/setup/samplesheet/samplesheet_{dataset}_{modality}.tsv` lists each
  library's FASTQs, one column per lane and read. Only datasets with samples
  to run get samplesheets.

## Sample names

A sample is named by the parts its library names share in the `sample_ids`
table (`25_13_1205_S25` and `25_13_1205_S41` give `25_13_1205`); sequencer
sample numbers (`S` + digits) are dropped. Add a `sample` column to the table
to name samples yourself. See `workflow/lib/rat_snake/samples.py`.

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `copy_config` | rule | Copy the run's config (after defaults are filled in) into the output folder. |
| `sample_status` | rule | Write every sample's libraries and status (run, existing or excluded, and why). |
| `make_samplesheet` | rule | List one dataset's FASTQs for one modality, one row per library with its lanes and reads. |
<!-- RULES:END -->
