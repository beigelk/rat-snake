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

# Parse command-line options
while getopts ":c:r:p:" opt; do
    case "$opt" in
        c)
            config="$OPTARG"
            ;;
        r)
        case "$OPTARG" in
            run|dry-run)
            runtype="$OPTARG"
            ;;
            *)
            log "ERROR" "Invalid runtype option: $OPTARG (allowed: run, dry-run)"
            exit 1
            ;;
        esac
        ;;
        p)
        case "$OPTARG" in
            local|slurm)
            profile="$OPTARG"
            ;;
            *)
            log "ERROR" "Invalid profile option: $OPTARG (allowed: local, slurm)"
            exit 1
            ;;
        esac
        ;;
        \?)
            log "ERROR" "Invalid option: -$OPTARG"
            echo "Usage: $0 -c <config.yaml> -r <dry-run|run>"
            exit 1
            ;;
        :)
            log "ERROR" "Option -$OPTARG requires an argument"
            exit 1
            ;;
    esac
done

# ───────────────────────────────────────────────
# Validate input
# ───────────────────────────────────────────────
config="${config:-}"
runtype="${runtype:-dry-run}"
# profile="${profile:-snakemake/profiles/local/}"

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
logdir="projects/${project_name}/${project}/logs"
mkdir -p ${logdir}

# ───────────────────────────────────────────────
# Make config dir
# ───────────────────────────────────────────────
configdir="projects/${project_name}/${project}/configs"
mkdir -p ${configdir}
cp ${config} ${configdir}

# ───────────────────────────────────────────────
# Run Snakemake
# ───────────────────────────────────────────────
log "INFO" "Running rat-snake pipeline"
log "INFO" "    Project     : $project_name ($project_version)"
log "INFO" "    Profile     : $profile"
log "INFO" "    Config file : $config"

dryrun_flag=""
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

profile_dir=""
case "$profile" in
    local)
        profile_dir="snakemake/profiles/local"
        ;;
    slurm)
        profile_dir="snakemake/profiles/slurm"
        ;;
    *)
        log "ERROR" "Unknown runtype: $runtype (allowed: local, slurm)"
        exit 1
        ;;
esac

# Main
log "INFO" "Starting rat-snake pipeline in '$runtype' mode"

script -q -e -c "
    snakemake \
        --snakefile snakemake/Snakefile \
        --profile '$profile_dir' \
        --configfile '$config' \
        $dryrun_flag \
        --rerun-triggers mtime \
        --printshellcmds \
        --use-singularity \
        --rerun-incomplete 
    " /dev/null | tee "${logdir}/rat-snake.log"

log "PASS" "sc-rat-snake Snakemake pipeline finished successfully."

# clean up
rm -rf .snakemake/tmp.*

    
