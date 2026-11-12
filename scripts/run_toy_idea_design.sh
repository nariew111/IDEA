#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${OUT_DIR:-${ROOT}/results/toy}"

if ! mkdir -p "${OUT}" 2>/dev/null; then
  OUT="${TMPDIR:-/tmp}/idea_toy_results"
  mkdir -p "${OUT}"
fi
PYTHON_BIN="${PYTHON_BIN:-python3}"

PYTHONPATH="${ROOT}/src" "${PYTHON_BIN}" -m idea.cli select \
  --species pig \
  --breeds toy \
  --num-bins 3 \
  --data-root "${ROOT}/examples/toy" \
  --output-dir "${OUT}"

echo "Selected SNPs:"
cat "${OUT}/toy_selected_snps.txt"
