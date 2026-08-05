#!/bin/bash

##################################################
################# USAGE EXAMPLE ##################
# bash rat-snake_run.sh -c config/project_configs.yaml -r run
##################################################
##################################################

set -euo pipefail
source snakemake/envs/conda/conda.sh

# ───────────────────────────────────────────────
# Logging helper
# ───────────────────────────────────────────────
log() {
  local level=$1
  local message=$2
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] [${level}] ${message}"
}

# ───────────────────────────────────────────────
# Parse input arguments
# ───────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)
            config="$2"
            shift 2
            ;;
        -r|--runtype)
            case "$2" in
                run|dry-run)
                    runtype="$2"
                    ;;
                *)
                    log "ERROR" "Invalid mode: $2 (allowed: run, dry-run)"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        -e|--executor)
            case "$2" in
                local|slurm)
                    executor="$2"
                    ;;
                *)
                    log "ERROR" "Invalid executor: $2 (allowed: local, slurm)"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            break
            ;;
    esac
done

# Everything remaining gets passed to snakemake
snakemake_args=("$@")

# ───────────────────────────────────────────────
# Validate input
# ───────────────────────────────────────────────
config="${config:-}"
runtype="${runtype:-dry-run}"
executor="${executor:-local}"

if [[ -z "${config:-}" ]]; then
    log "ERROR" "Missing -c argument"
    echo "Usage: $0 -c <path_to_config.yaml>"
    exit 1
fi

if [[ ! -f "$config" ]]; then
    log "ERROR" "File not found: $config"
    exit 1
fi

# ───────────────────────────────────────────────
# Check for yq
# ───────────────────────────────────────────────
if ! command -v yq >/dev/null 2>&1; then
    log "ERROR" "'yq' is required but not installed."
    exit 1
fi

# ───────────────────────────────────────────────
# Parse project name from JSON
# ───────────────────────────────────────────────
log "INFO" "Reading project configuration from: $config"

# Fail if JSON is invalid
if ! yq -e '.' "$config" >/dev/null 2>&1; then
    log "ERROR" "Invalid JSON file: $config"
    exit 1
fi

project_name=$(yq -r '.project_info.name // empty' "$config" | xargs)
project_version=$(yq -r '.project_info.version // "0.1.0"' "$config" | xargs)

if [[ -z "$project_name" ]]; then
    log "ERROR" "'project_name' not found or is empty in $config"
    exit 1
fi

project="${project_name}_${project_version}"

# ───────────────────────────────────────────────
# Make log dir
# ───────────────────────────────────────────────
logdir="projects/${project_name}/${project}/.logs"
mkdir -p "${logdir}"

# ───────────────────────────────────────────────
# Make config dir
# ───────────────────────────────────────────────
configdir="projects/${project_name}/${project}/configs"
mkdir -p "${configdir}/input"
cp "${config}" "${configdir}/input"

# ───────────────────────────────────────────────
# Run Snakemake
# ───────────────────────────────────────────────
log "INFO" "Running rat-snake pipeline"
log "INFO" "    Project      : $project_name ($project_version)"
log "INFO" "    Executor     : $executor"
log "INFO" "    Config file  : $config"

case "$runtype" in
    run)
        dryrun_flag=""
        ;;
    dry-run)
        dryrun_flag="--dry-run"
        ;;
    *)
        log "ERROR" "Unknown runtype: $runtype (allowed: run, dry-run)"
        exit 1
        ;;
esac

case "$executor" in
    local)
        profile_dir="snakemake/profiles/local"
        ;;
    slurm)
        profile_dir="snakemake/profiles/slurm"
        ;;
    *)
        log "ERROR" "Unknown executor: $executor (allowed: local, slurm)"
        exit 1
        ;;
esac

PROJECT_ROOT=$(pwd)
export R_ENVIRON_USER=${PROJECT_ROOT}/.Renviron
export R_PROFILE_USER=${PROJECT_ROOT}/.Rprofile

export FLEXIBLAS=NETLIB

# ───────────────────────────────────────────────
# Build and run snakemake via temp script
# ───────────────────────────────────────────────
log "INFO" "Starting rat-snake pipeline in '$runtype' mode"
log "INFO" "Passing the following additional flags to Snakemake: ${snakemake_args[@]}"

tmp_script=$(mktemp /tmp/rat-snake-XXXXXX.sh)
trap 'rm -f "$tmp_script"' EXIT

cat > "$tmp_script" << EOF
#!/bin/bash
snakemake \\
    --snakefile snakemake/Snakefile \\
    --profile '$profile_dir' \\
    --configfile '$config' \\
    ${dryrun_flag:+--dry-run} \\
    --printshellcmds \\
    --use-envmodules \\
    ${snakemake_args[@]+"${snakemake_args[@]}"}
EOF

script -q -e -c "bash '$tmp_script'" /dev/null | tee "${logdir}/rat-snake.log"

log "PASS" "rat-snake Snakemake pipeline finished successfully."

# clean up
rm -rf .snakemake/tmp.*
