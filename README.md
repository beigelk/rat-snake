# rat-snake

**R**NA, **A**DT and **T**CR snake: a Snakemake pipeline for 10x Genomics
CITE-seq + TCR single-cell data, from FASTQs (or finished Cell Ranger runs)
to clustered, annotated populations, marker tables, plots and a report.

It is built on [snake-core](../snake-core) (in `core/`), the shared base of
the snake pipeline suite (sunbeam-snake uses the same core): the controller,
logging, output paths, software modes, R script set-up and the report
builder.

## Modules

They always run in this order; each needs the one before it.

| Module | What it does |
|---|---|
| [setup](workflow/modules/setup/README.md) | Copies the config, gives every sample a status (run, existing, excluded), lists FASTQs |
| [cellranger](workflow/modules/cellranger/README.md) | Cell Ranger multi per sample (or finished runs), metrics, web summaries |
| [preprocess](workflow/modules/preprocess/README.md) | SoupX, DoubletFinder, dsb, scRepertoire clonotypes, filtered RNA + ADT object per sample |
| [integrate](workflow/modules/integrate/README.md) | Merge, Harmony on RNA and ADT, WNN and UMAP |
| [annotate](workflow/modules/annotate/README.md) | Sample metadata and clonotypes on the cells; per population: clusters, subsets, markers, clustree |
| [visualization](workflow/modules/visualization/README.md) | Metadata UMAPs and QC plots; dot plots and marker plots per marker set |
| [report](workflow/modules/report/README.md) | One report of figures and tables as HTML, Word and Notion-ready Word |

`modules` in the config can stop the run early (e.g. `[setup, cellranger]`).

## Getting the code

The pipeline needs its `core/` submodule.

### From GitHub

```bash
git clone --recurse-submodules -b restructure https://github.com/<you>/rat-snake.git
```

### From the git bundles

`rat-snake.bundle` and `snake-core.bundle` are git repositories in one file
each, history included. In your existing rat-snake clone:

```bash
git fetch /path/to/rat-snake.bundle restructure:restructure
git checkout restructure
git config submodule.core.url /path/to/snake-core.bundle
git -c protocol.file.allow=always submodule update --init   # git blocks local submodules otherwise
```

To push: `git push -u origin restructure`. `.gitmodules` uses the relative
URL `../snake-core.git`, so on GitHub the submodule resolves to the
`snake-core` repo next to `rat-snake`; push snake-core there first
(`cd core && git push <snake-core remote> main`), then
`git config --unset submodule.core.url`.

## Set up (once)

```bash
bash envs/create_module_envs.sh       # one conda env per module env file (rat-preprocess-r, ...)
```

The env that runs Snakemake itself is created by the controller the first
time it runs (`envs/conda_setup.sh`). Cell Ranger is not a conda package:
set `cellranger.path` or `software.envmodules.cellranger` (see
[cellranger](workflow/modules/cellranger/README.md)).

## Run

```bash
cp config/config.yaml config/my_project.yaml       # then edit it
bash rat-snake.sh -c config/my_project.yaml               # dry run (default)
bash rat-snake.sh -c config/my_project.yaml -r run        # run here
bash rat-snake.sh -c config/my_project.yaml -r run -e slurm
bash rat-snake.sh -c config/my_project.yaml -r run -- --until harmony   # anything after -- goes to snakemake
bash rat-snake.sh -c config/my_project.yaml -r jobs       # this project's Slurm jobs, with rule and sample
```

### Slurm job names

The Slurm executor gives every job of a run the same name and uses that name
to track the jobs, so the name cannot differ per rule or sample. With
`-e slurm` the name starts with the pipeline and project
(`rat-snake_<name>_<version>_<run id>`), and the rule and wildcards (dataset,
sample, population) go in the job's comment, e.g.
`rule_cellranger_multi_wildcards_d1_25_13_1205`. `-r jobs` lists the
project's queued and running jobs with that comment as a "RULE" column. With
plain `squeue`, add the comment column (`%k`):

```bash
squeue --me --format "%.10i %.9T %.10M %.40j %k"
```

Each job's Slurm log is under `.snakemake/slurm_logs/rule_<rule>/<wildcards>/<jobid>.log`.

Every option is in `config/config.yaml` and described, with its default, in
`config/config.schema.yaml`; the config is checked against the schema at
start-up. Each dataset needs a `sample_ids` table (tab-separated, columns
`gex`, `csp`, `vdj` with each library's FASTQ sample name, and an optional
`sample` column) and its FASTQ folders.

Outputs go to `output/{name}_{version}/`:

```
config/config.yaml       the config with defaults filled in
results/{module}/...     every table, object and figure
reports/                 report.html, report.docx, report_notion.docx
logs/{module}/           one log per job (logs/latest.log is the newest controller log)
benchmarks/{module}/     time and memory per job
debug/{module}/          saved Snakemake objects (debug: true)
```

## Disk space

Nothing the pipeline writes goes to your home folder: results, logs and
Snakemake's state are in `output/` and `.snakemake/`, temporary files in
`tmp/`, caches (Snakemake's source cache, pip, the Apptainer image cache) in
`.snakemake/cache/` (set `CACHE_DIR_OVERRIDE` to use another folder), and
Cell Ranger runs in `cellranger.scratch_dir`. So keep the pipeline folder on
a project or scratch file system, not in home (on our cluster:
`/scr1/users/$USER`, 20 TB). To keep the outputs there while the code stays
elsewhere, set `project.output_root` to a folder under it (write the path
out; `$USER` is not expanded in the config).

Conda decides where envs and downloaded packages go (often `~/.conda`).
The module envs take several GB, so on a small or shared home, point conda
elsewhere before running `envs/create_module_envs.sh`, in `~/.condarc`:

```yaml
envs_dirs: [/scr1/users/<username>/conda/envs]
pkgs_dirs: [/scr1/users/<username>/conda/pkgs]
```

The controller warns when conda still uses home.

## Debugging

**Skip Cell Ranger.** Point `cellranger.existing_outs` at finished Cell
Ranger multi runs and those samples are never re-run:

```yaml
cellranger:
  existing_outs:
    root: /data/old_run/cellranger              # {root}/{dataset}/{sample} or {root}/{sample}
    samples:
      25_13_1205: /data/runs/25_13_1205/outs    # one sample, by sample or dataset/sample
```

**Re-run an R step by hand.** With `debug: true`, every R job saves its
Snakemake object to `debug/{module}/{job}.rds`. In R, from the repo root:

```r
Sys.setenv(SNAKE_CORE = "core")
snakemake <- readRDS("output/my_project_v0/debug/annotate/subset__lymphocytes.rds")
source("workflow/modules/annotate/scripts/subset.R")   # or step through it
```

## Software environments

`software.mode` picks where each rule's software comes from:

| Mode | Software |
|---|---|
| `conda` (default) | Named envs from `envs/create_module_envs.sh`, e.g. `rat-annotate-r` |
| `conda_yaml` | Snakemake builds each env from its env file |
| `container` | Singularity/Apptainer images `{registry}/rat-{module}[-{variant}]:{image_tag}`, built by GitHub Actions (`.github/workflows/images.yml`) |
| `none` | Whatever is on PATH |

Each module has its env files next to its rules: `env_r.yaml` for R steps
(preprocess also `env_doublet.yaml`, report `env.yaml`). R packages conda
does not have go in `<env>.post.R` (visualization installs scplotter and
plotthis from CRAN this way). `software.r_library` puts an extra R library
(e.g. a renv library) first on every script's library path.

## Tests

```bash
bash tests/run_tests.sh       # unit tests, rule docs check, dry runs (no R needed)
bash tests/run_synthetic.sh   # every module on a small synthetic project (conda envs)
```

The synthetic project (`tests/make_synthetic.py`) has three samples in two
datasets, one more sample that is excluded for lacking a VDJ library, four
cell types, ADT and TCR data, as finished Cell Ranger runs.
`tests/test_synthetic_outputs.py` checks the results.

## Layout

```
rat-snake.sh                 controller (core/shell/controller.sh)
config/config.yaml           example config; config.schema.yaml documents every option
profiles/{local,slurm}/      Snakemake profiles (resources per rule)
envs/                        Snakemake env and module env set-up
workflow/Snakefile           reads the config and includes the modules
workflow/lib/rat_snake/      samples, populations and their paths (Python)
workflow/lib/R/              R helpers shared by scripts
workflow/modules/{module}/   {module}.smk, targets.py, scripts/, env files, README.md
core/                        snake-core (submodule)
tests/                       synthetic data, unit tests, test scripts
```

## Documentation

Each module's README describes its outputs and options and ends with a rule
table generated from the rule docstrings
(`python core/bin/update_rule_docs.py workflow/modules`). Every script starts
with a header saying what it does, and every function is documented.
[MIGRATION.md](MIGRATION.md) maps the old layout and config keys to the new
ones.
