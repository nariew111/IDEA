#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: bash configs/idea_gp.example.sh PLINK_PREFIX PHENOTYPE_FILE OUTPUT_PREFIX" >&2
  exit 2
fi

# This is an example invocation, not the manuscript's locked benchmark settings.
# Install the external tools required by the chosen model before running.
idea-gp run \
  --bfile "$1" \
  --pheno "$2" \
  --out "$3" \
  --folds 5 \
  --inner_cv_folds 5 \
  --trials 20 \
  --method rrblup \
  --optuna_seed 2026
