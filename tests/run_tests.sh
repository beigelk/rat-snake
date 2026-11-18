#!/usr/bin/env bash
# Tests that need no R and no internet: unit tests, the rule-docs check and
# dry runs of both synthetic configs (finished Cell Ranger runs, and running
# Cell Ranger).
#
# Usage: bash tests/run_tests.sh   (from anywhere, with snakemake, pandas,
#        pyyaml and pytest on PATH)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

python tests/make_synthetic.py --out tests/data
python -m pytest -q tests/test_lib.py core/tests
python core/bin/update_rule_docs.py workflow/modules --check
for cfg in config_existing config_run; do
    snakemake -s workflow/Snakefile --configfile "tests/data/${cfg}.yaml" -n --cores 1 -q
done
snakemake -s workflow/Snakefile --configfile tests/data/config_existing.yaml -n --cores 1 -q \
    --config modalities='[gex,csp]'
echo "All tests passed"
