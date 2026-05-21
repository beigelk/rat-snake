#!/usr/bin/env bash

RUN_ENV="snakemake-v9.17.2"

SMK_V="9.17.2"

REQ_PKGS=(
    "snakemake=${SMK_V}"
    "snakemake-executor-plugin-slurm"
    "yq"
)

# Function for checking if a conda env exists
# ---------------------------------------------------------------------------------------
# Source - https://stackoverflow.com/a/70598193
# Posted by kevinnlsamuel, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-23, License - CC BY-SA 4.0

# Usage
# if find_in_conda_env ".*RUN_ENV.*" ; then
#    conda activate RUN_ENV
# else 

find_in_conda_env() {
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}
# ---------------------------------------------------------------------------------------

# Initialize conda for this shell
source "$(conda info --base)/etc/profile.d/conda.sh"

if find_in_conda_env "${RUN_ENV}" ; then
   conda activate "${RUN_ENV}"
else

    # Ensure channels are in the right order
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --set channel_priority strict

    conda create -n "${RUN_ENV}" \
        -c conda-forge -c bioconda \
        "${REQ_PKGS[@]}" \
        -y

    conda activate "${RUN_ENV}"
fi



# Core paths
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$CONDA_PREFIX/lib:${LIBRARY_PATH:-}"
export CPATH="$CONDA_PREFIX/include:${CPATH:-}"

# FreeType
export CPATH="$CONDA_PREFIX/include/freetype2:${CPATH:-}"

# Java (rJava requirement)
export JAVA_HOME="$CONDA_PREFIX"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/jvm/lib/server:$CONDA_PREFIX/lib/jvm/lib:${LD_LIBRARY_PATH:-}"
export CPATH="$CONDA_PREFIX/lib/jvm/include:$CONDA_PREFIX/lib/jvm/include/linux:${CPATH:-}"