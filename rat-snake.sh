#!/usr/bin/env bash

# Snakemake pipeline for 10x CITE-seq + TCR single-cell data
# RNA (r), ADT (a) and TCR (t) snake: Cell Ranger multi, ambient RNA, doublets,
# ADT normalisation, clonotypes, Harmony/WNN integration, clustering, markers,
# plots and a report.

##################################################
################# USAGE EXAMPLE ##################
# bash rat-snake.sh \
#   -c config/config.yaml \
#   -r run \
#   -e slurm
#
# Anything after -- goes straight to snakemake:
# bash rat-snake.sh -c config/config.yaml -r run -- --until harmony
##################################################

PIPELINE_NAME="rat-snake"
SNAKEFILE="workflow/Snakefile"
PROFILES_DIR="profiles"

cd "$(dirname "${BASH_SOURCE[0]}")"

if [[ ! -f core/shell/controller.sh ]]; then
    echo "core/ is empty. Run: git submodule update --init" >&2
    exit 1
fi

source core/shell/controller.sh
snake_main "$@"
