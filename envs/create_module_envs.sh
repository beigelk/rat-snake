#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# create_module_envs.sh: build (or update) the named conda envs of each module.
#
# For software.mode: conda. Each env file of a module (env.yaml, env_r.yaml,
# env_doublet.yaml, ...) becomes the env named in its "name:" field (e.g.
# rat-preprocess-r). Rules then use the env by name, so Snakemake never
# rebuilds it, and you can work in it directly:
#
#   conda activate rat-annotate-r
#
# When an env file has a post-install R script next to it (env_r.yaml ->
# env_r.post.R), it runs in the env afterwards, for R packages conda does not
# have (scplotter and plotthis in the visualization module).
#
# Usage:
#   bash envs/create_module_envs.sh                         # every module, create or update
#   bash envs/create_module_envs.sh preprocess annotate     # just these modules
#   bash envs/create_module_envs.sh --missing preprocess    # only create envs that do not exist yet
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source core/shell/log.sh
source "$(conda info --base)/etc/profile.d/conda.sh"

solver=conda
command -v mamba >/dev/null 2>&1 && solver=mamba

missing_only=false
modules=()
for arg in "$@"; do
    case "$arg" in
        --missing) missing_only=true ;;
        *) modules+=("$arg") ;;
    esac
done
if [[ ${#modules[@]} -eq 0 ]]; then
    for d in workflow/modules/*/; do
        compgen -G "${d}env*.yaml" >/dev/null && modules+=("$(basename "$d")")
    done
fi

# env_name FILE: the "name:" field of an env file (no yq needed)
env_name() {
    awk '/^name:/ {print $2; exit}' "$1"
}

for m in "${modules[@]}"; do
    compgen -G "workflow/modules/${m}/env*.yaml" >/dev/null || die "No env files for module '${m}'"
    for env_file in workflow/modules/"${m}"/env*.yaml; do
        name=$(env_name "$env_file")
        post="${env_file%.yaml}.post.R"
        if conda env list | awk '{print $1}' | grep -qx "$name"; then
            if [[ "$missing_only" == true ]]; then
                log "INFO" "${name} exists; skipping"
                continue
            fi
            log "INFO" "Updating ${name} from ${env_file}"
            "$solver" env update -n "$name" -f "$env_file" --prune
        else
            log "INFO" "Creating ${name} from ${env_file}"
            "$solver" env create -n "$name" -f "$env_file"
        fi
        if [[ -f "$post" ]]; then
            log "INFO" "Running ${post} in ${name}"
            conda run -n "$name" Rscript "$post"
        fi
    done
done
log "PASS" "Module envs ready: ${modules[*]}"
