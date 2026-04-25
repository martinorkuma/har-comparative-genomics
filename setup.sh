#!/usr/bin/env bash
# One-time setup. Idempotent: safe to re-run.
#
# Usage:
#   bash setup.sh
#
# What it does:
#   1. Creates the conda env `har-ml` if it doesn't exist.
#   2. Verifies key CLI tools resolve.
#   3. Creates the data/, outputs/, logs/ directories if missing.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
cd "$ROOT"

# ---- conda env ----
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install Miniconda or Anaconda first." >&2
    exit 1
fi

# Source conda's shell hook so `conda activate` works inside the script
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx 'har-ml'; then
    echo "[setup] conda env 'har-ml' already exists; skipping creation."
else
    echo "[setup] creating conda env 'har-ml' from environment.yml ..."
    conda env create -f environment.yml
fi

conda activate har-ml

# ---- tool check ----
echo "[setup] checking required tools ..."
for tool in bedtools liftOver bigBedToBed bigWigAverageOverBed python; do
    if command -v "$tool" &> /dev/null; then
        echo "  ok: $tool ($($tool --version 2>&1 | head -n1))"
    else
        echo "  MISSING: $tool" >&2
        exit 1
    fi
done

# ---- python imports check ----
python - <<'PY'
import importlib, sys
mods = ["numpy", "pandas", "scipy", "openpyxl", "sklearn", "matplotlib", "seaborn",
        "yaml", "pybedtools", "pyfaidx", "pyBigWig", "shap"]
missing = []
for m in mods:
    try:
        importlib.import_module(m)
    except ImportError:
        missing.append(m)
if missing:
    print("MISSING python modules:", missing, file=sys.stderr)
    sys.exit(1)
print("[setup] all python modules import cleanly.")
PY

# ---- directories ----
mkdir -p data/raw data/processed outputs/figures outputs/tables outputs/models logs

echo
echo "[setup] done. activate the env with:  conda activate har-ml"
echo "[setup] then run:                    bash scripts/run_all.sh"
