#!/usr/bin/env bash
# The env that runs Snakemake itself (not the analysis software; each module
# has its own env files). Created once, then just activated.

RUN_ENV="smk-v9.17.2__rat__v1"

SMK_V="9.17.2"

REQ_PKGS=(
    "snakemake=${SMK_V}"
    "snakemake-executor-plugin-slurm"
    "python=3.11"
    "pyyaml"
    "jsonschema"
    "pandas"
    "yq"
)

source core/shell/conda_setup.sh
