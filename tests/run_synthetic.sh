#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# run_synthetic.sh: run the whole pipeline on a small synthetic project.
#
# Three samples in two datasets (plus one excluded for lacking a VDJ library),
# four cell types, ADT and TCR data. Cell Ranger is not run: the samples are
# finished Cell Ranger runs found through cellranger.existing_outs.
#
# Steps: check that core/ matches this commit; activate (or create) the env
# that runs Snakemake; create any missing module envs; write the synthetic
# data; run every module; check the results with pytest.
#
# Usage (from anywhere):
#   bash tests/run_synthetic.sh              # conda envs (software.mode: conda)
#   bash tests/run_synthetic.sh -j 8         # cores for Snakemake
#   bash tests/run_synthetic.sh -e slurm     # run the jobs through Slurm
#   bash tests/run_synthetic.sh -m none      # use the active env for everything
#   bash tests/run_synthetic.sh --keep       # keep earlier output (resume)
#
# Results go to tests/data/output/synthetic_v0/ (report in reports/).
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source core/shell/log.sh 2>/dev/null || {
    echo "core/ is empty. Run: git -c protocol.file.allow=always submodule update --init" >&2
    exit 1
}

mode="conda" cores="" executor="local" keep=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mode)     mode="$2"; shift 2 ;;
        -j|--cores)    cores="$2"; shift 2 ;;
        -e|--executor) executor="$2"; shift 2 ;;
        --keep)        keep=true; shift ;;
        -h|--help)     sed -n '2,21p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)             die "Unknown option: $1 (see --help)" ;;
    esac
done
case "$mode" in conda|none) ;; *) die "-m must be conda or none" ;; esac

# ── 1. core/ must be the snake-core commit this rat-snake expects ────────────
if git rev-parse --git-dir >/dev/null 2>&1; then
    status=$(git submodule status core 2>/dev/null || true)
    case "${status:0:1}" in
        -) die "core/ is not checked out. Run: git -c protocol.file.allow=always submodule update --init" ;;
        +) die "core/ is on $(git -C core rev-parse --short HEAD) but this commit needs $(git ls-tree HEAD core | awk '{print substr($3,1,7)}').
       Update snake-core first, then run:
         git submodule sync && git -c protocol.file.allow=always submodule update --init" ;;
    esac
fi
log "PASS" "core/ matches this commit ($(git -C core log --oneline -1 2>/dev/null || echo unknown))"

# ── 2. Software ──────────────────────────────────────────────────────────────
controller_flags=()
if [[ "$mode" == conda ]]; then
    command -v conda >/dev/null 2>&1 || die "conda is not on PATH (load your conda module, or use -m none)"
    log "INFO" "Activating the Snakemake env (created the first time; this can take a while)"
    source envs/conda_setup.sh
    command -v pytest >/dev/null 2>&1 || conda install -y -c conda-forge pytest
    log "INFO" "Creating missing module envs"
    bash envs/create_module_envs.sh --missing
else
    controller_flags=(-s)
fi

# ── 3. Data and config ───────────────────────────────────────────────────────
python tests/make_synthetic.py --out tests/data
config=tests/data/config_synthetic.yaml
sed -E "s/^([[:space:]]*mode:[[:space:]]*)none/\1${mode}/" tests/data/config_existing.yaml > "$config"
[[ "$keep" == true ]] || rm -rf tests/data/output/synthetic_v0

# ── 4. Run and check ─────────────────────────────────────────────────────────
snakemake_args=()
[[ -n "$cores" ]] && snakemake_args+=(--cores "$cores")
log "START" "Running every module on ${config}"
bash rat-snake.sh -c "$config" -r run -e "$executor" "${controller_flags[@]}" -- "${snakemake_args[@]}" \
    || die "Pipeline failed (logs in tests/data/output/synthetic_v0/logs/)"
python -m pytest -q tests/test_synthetic_outputs.py
log "PASS" "Synthetic run passed; report: tests/data/output/synthetic_v0/reports/report.html"
