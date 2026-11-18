#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# run_cellranger_multi.sh: run Cell Ranger multi for one sample.
#
# Called by rule cellranger_multi. With a scratch folder (cellranger.scratch_dir)
# Cell Ranger runs there and the finished run is copied to OUTDIR, as the old
# rule did on the cluster; without one it runs next to OUTDIR. The large
# SC_MULTI_CS folder is deleted unless KEEP_CS is true.
#
# Usage:
#   run_cellranger_multi.sh CONFIG_CSV SAMPLE OUTDIR THREADS MEM_GB SCRATCH KEEP_CS CELLRANGER
#     SCRATCH     folder or "" (environment variables such as $USER are expanded)
#     CELLRANGER  cellranger executable ("" = cellranger on PATH)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

csv=$(realpath "$1"); sample=$2; outdir=$(realpath -m "$3")
threads=$4; mem_gb=$5; scratch=$6; keep_cs=$7; cellranger=${8:-cellranger}
[[ -n "$cellranger" ]] || cellranger=cellranger

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$1] $2"; }

command -v "$cellranger" >/dev/null 2>&1 || {
    log ERROR "'$cellranger' not found. Set cellranger.path or software.envmodules.cellranger in the config."
    exit 1
}
log INFO "Cell Ranger: $("$cellranger" --version 2>&1 | head -1)"

if [[ -n "$scratch" ]]; then
    scratch=$(envsubst <<< "$scratch" 2>/dev/null || eval "echo \"$scratch\"")
    workdir="${scratch%/}/rat-snake_${SLURM_JOB_ID:-$$}_${sample}"
else
    workdir=$(dirname "$outdir")
fi
mkdir -p "$workdir"
log INFO "Working folder: $workdir"

cd "$workdir"
rm -rf "${workdir:?}/${sample}"   # a stale pipestance from an earlier failed run
"$cellranger" multi \
    --id="$sample" \
    --csv="$csv" \
    --localcores="$threads" \
    --localmem="$mem_gb" \
    --disable-ui

if [[ "$keep_cs" != "true" ]]; then
    rm -rf "${workdir:?}/${sample}/SC_MULTI_CS"
    log INFO "Deleted SC_MULTI_CS"
fi

if [[ -n "$scratch" ]]; then
    mkdir -p "$(dirname "$outdir")"
    rm -rf "$outdir"
    cp -r "${workdir}/${sample}" "$outdir"
    log INFO "Copied results to $outdir"
    rm -rf "$workdir"
    log INFO "Removed $workdir"
fi
log DONE "cellranger multi finished for ${sample}"
