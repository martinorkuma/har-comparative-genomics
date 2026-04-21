#!/usr/bin/env bash
"""
setup_env.sh
One-time environment setup for the HAR comparative genomics project.
Builds a conda env containing both the bioinformatics binaries (bedtools,
liftOver, UCSC tools) and all Python dependencies, then creates the
project directory tree.

Works on WSL (Ubuntu), Linux, and macOS.

Usage:
   bash setup_env.sh              # full install
   bash setup_env.sh --minimal    # skip heavy extras (shap, jupyterlab, reportlab)
                                  # useful for CI or Phase 1-only runs
   bash setup_env.sh --update     # update an existing env in place

 Author: Martin Orkuma
 Project: HAR Comparative Genomics

 Debugging and code assistance were provided by ChatGPT (GPT 5.4 Thinking) and Claude (Opus 4.7).
"""

set -euo pipefail

# ---- Parse flags ----
MODE="full"
for arg in "$@"; do
    case "$arg" in
        --minimal) MODE="minimal" ;;
        --update)  MODE="update" ;;
        -h|--help)
            sed -n '2,20p' "$0"; exit 0 ;;
        *)
            echo "Unknown flag: $arg"; echo "Run '$0 --help' for usage."; exit 1 ;;
    esac
done

# ---- Resolve repo root so the script works no matter where it's invoked ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
cd "$REPO_ROOT"

echo "=================================================="
echo " HAR Comparative Genomics — Environment Setup"
echo " Repo:  $REPO_ROOT"
echo " Mode:  $MODE"
echo "=================================================="

ENV_NAME="har-ml"
ENV_YML="environment/environment.yml"
REQ_TXT="environment/requirements.txt"

# ---- Check conda is installed ----
# We need conda (or mamba) because bedtools, liftOver, and bigWigAverageOverBed
# are compiled binaries distributed via bioconda — pip can't install them.
if ! command -v conda >/dev/null 2>&1; then
    cat <<EOF
ERROR: conda not found.

This project depends on bedtools, liftOver, and UCSC utilities that are
only distributed through bioconda. Install Miniconda first:

  https://docs.conda.io/en/latest/miniconda.html

On WSL/Ubuntu the quick path is:

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh

then restart your shell and re-run this script.
EOF
    exit 1
fi

# Prefer mamba if available — it resolves bioconda deps 5-10x faster than conda.
SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then
    SOLVER="mamba"
    echo "  [info] mamba detected — using it for faster solving"
fi

# ---- Verify declarative env files exist ----
for f in "$ENV_YML" "$REQ_TXT"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found. Did you run this from the project root?"
        exit 1
    fi
done

# ---- Create, update, or warn about existing env ----
# conda env list output has the env name as the first whitespace token per line.
env_exists () {
    conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"
}

if [ "$MODE" = "update" ] || env_exists; then
    if env_exists; then
        echo "  [info] env '$ENV_NAME' exists — updating in place"
        $SOLVER env update -n "$ENV_NAME" -f "$ENV_YML" --prune
    else
        echo "ERROR: --update requested but env '$ENV_NAME' does not exist."
        exit 1
    fi
else
    echo "  [info] creating env '$ENV_NAME' from $ENV_YML"
    $SOLVER env create -f "$ENV_YML"
fi

# ---- If --minimal, uninstall the heavy extras we don't need for Phase 1-2 ----
# We still install them by default because Phase 3-4 need them; this flag is
# an escape hatch for CI or for getting a quick working env on a slow laptop.
if [ "$MODE" = "minimal" ]; then
    echo "  [info] --minimal: removing shap, jupyterlab, reportlab, pypdf"
    # 'conda run' executes inside the env without requiring shell activation.
    conda run -n "$ENV_NAME" pip uninstall -y shap jupyterlab ipykernel reportlab pypdf || true
fi

# ---- Create project directory tree ----
# These match the paths referenced by Phase 1 scripts (data/raw, data/processed,
# logs/) and anticipate Phase 2-4 needs (outputs/figures etc.).
# .gitkeep files let git track the empty dirs so collaborators don't have to
# recreate them before running the pipeline.
echo "  [info] creating project directories"
mkdir -p data/raw data/processed
mkdir -p logs
mkdir -p outputs/figures outputs/tables outputs/models outputs/reports

touch data/raw/.gitkeep data/processed/.gitkeep
touch logs/.gitkeep
touch outputs/figures/.gitkeep outputs/tables/.gitkeep
touch outputs/models/.gitkeep outputs/reports/.gitkeep

# ---- Quick smoke test: confirm the key binaries and packages are reachable ----
echo
echo "  [info] smoke test — verifying key tools inside '$ENV_NAME'"
conda run -n "$ENV_NAME" bash -c '
    set -e
    for tool in bedtools liftOver; do
        if command -v "$tool" >/dev/null 2>&1; then
            printf "    [OK]  %-25s %s\n" "$tool" "$(command -v "$tool")"
        else
            printf "    [MISS] %s\n" "$tool"; exit 1
        fi
    done
    python - <<PY
import importlib
for p in ["pandas", "numpy", "sklearn"]:
    m = importlib.import_module(p)
    print(f"    [OK]  {p:25s} {getattr(m, \"__version__\", \"?\")}")
PY
'

# ---- Register Jupyter kernel so the env shows up in notebooks (full mode only) ----
if [ "$MODE" = "full" ]; then
    echo "  [info] registering Jupyter kernel 'har-ml'"
    conda run -n "$ENV_NAME" python -m ipykernel install \
        --user --name "$ENV_NAME" --display-name "Python (har-ml)" >/dev/null 2>&1 || \
        echo "    [warn] kernel install skipped (non-fatal)"
fi

cat <<EOF

==================================================
 Setup complete!
==================================================

Next steps:

  1) Activate the environment:
       conda activate $ENV_NAME

  2) Run Phase 1 (data acquisition):
       bash scripts/01_data_acquisition/run_phase1.sh

  3) When you add a new Python dep, edit environment/requirements.txt
     and run:
       bash setup_env.sh --update

EOF
