#!/usr/bin/env bash
# run_all.sh — run the full pipeline end-to-end.
#
# Idempotent: each stage skips work whose outputs already exist.
# To force a stage to re-run, delete its outputs first.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &> /dev/null && pwd)"
cd "$ROOT"

# ---- env check ----
if [[ "${CONDA_DEFAULT_ENV:-}" != "har-ml" ]]; then
    echo "WARN: conda env 'har-ml' is not active. Run:  conda activate har-ml" >&2
fi

mkdir -p logs

stage() {
    local label="$1"; shift
    local logfile="logs/${label}.log"
    echo
    echo "================================================================="
    echo "  $label"
    echo "================================================================="
    echo "  log: $logfile"
    "$@" 2>&1 | tee "$logfile"
}

stage "01_acquire"       bash   scripts/01_acquire.sh
stage "02_build_features" python scripts/02_build_features.py
stage "03_classify"       python scripts/03_classify.py
stage "04_interpret"      python scripts/04_interpret.py

echo
echo "================================================================="
echo "  pipeline complete."
echo "================================================================="
echo "  poster figures: outputs/figures/"
echo "  tables:         outputs/tables/"
echo "  models:         outputs/models/"
