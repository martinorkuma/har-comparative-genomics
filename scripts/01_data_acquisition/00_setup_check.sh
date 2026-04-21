#!/usr/bin/env bash
"""
00_setup_check.sh
Verifies required command-line tools and Python packages are available before the Phase 1 pipeline runs. 
Fails fast with actionable error messages rather than halfway through a multi-gigabyte download.

If the user isn't in the 'har-ml' conda env, still proceed, and only fail if an actual tool is missing.
"""
set -euo pipefail

echo "=== Phase 1 setup check ==="
echo

# ---- 1. Nudge toward the canonical conda env ----
# CONDA_DEFAULT_ENV is set by `conda activate`. If it's not 'har-ml' we warn
# but continue, because users legitimately might run this on a shared cluster,
# in a docker image, or with apt-installed bedtools.
if [ "${CONDA_DEFAULT_ENV:-}" != "har-ml" ]; then
    cat <<EOF
NOTE: the 'har-ml' conda environment is not active.
      (current CONDA_DEFAULT_ENV: '${CONDA_DEFAULT_ENV:-<none>}')

The reproducible way to run this project is:
    bash setup_env.sh          # one-time, from the project root
    conda activate har-ml      # before each session

Proceeding anyway — if the tools below are on your PATH from another
install route (apt, homebrew, a different conda env) we'll use them.

EOF
fi

MISSING=()

check_tool () {
    local tool="$1"
    local install_hint="$2"
    if command -v "$tool" >/dev/null 2>&1; then
        printf "  [OK]   %-22s %s\n" "$tool" "$(command -v "$tool")"
    else
        printf "  [MISS] %-22s install: %s\n" "$tool" "$install_hint"
        MISSING+=("$tool")
    fi
}

# ---- 2. Check command-line tools ----
# Install hints prefer conda (the canonical path) but mention the apt/binary
# fallback where one exists, so users on locked-down systems aren't stuck.
echo "Command-line tools:"
check_tool curl                 "apt install curl"
check_tool gunzip               "apt install gzip"
check_tool bedtools             "conda install -c bioconda bedtools   OR   apt install bedtools"
check_tool liftOver             "conda install -c bioconda ucsc-liftover   OR   grab from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
check_tool bigWigAverageOverBed "conda install -c bioconda ucsc-bigwigaverageoverbed   (needed in Phase 2)"

# ---- 3. Check Python packages ----
# We only check the two we strictly need in Phase 1. Phase 2+ deps
# (scikit-learn, shap, biopython, ...) are verified by each phase's
# own setup check, so a minimal install can still clear Phase 1.
echo
echo "Python packages:"
python3 - <<'PY'
import importlib, sys
pkgs = ["pandas", "numpy"]
missing = []
for p in pkgs:
    try:
        m = importlib.import_module(p)
        print(f"    [OK]  {p:25s} {getattr(m, '__version__', '?')}")
    except ImportError:
        print(f"  [MISS] {p:10s} install: pip install {p}   (or: bash setup_env.sh --update)")
        missing.append(p)
sys.exit(1 if missing else 0)
PY

# ---- 4. Summarise and exit ----
if [ ${#MISSING[@]} -gt 0 ]; then
    cat <<EOF

Missing command-line tools: ${MISSING[*]}

Fastest fix — from the project root, run:
    bash setup_env.sh                 # first-time setup of 'har-ml'
    bash setup_env.sh --update        # env exists but has drifted

Then re-run this check. If you installed tools via apt or manually,
make sure they're on your PATH and re-run.
EOF
    exit 1
fi

echo
echo "All required tools found. Ready to proceed."