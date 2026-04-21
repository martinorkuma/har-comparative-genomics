#!/usr/bin/env bash
# 00_setup_check.sh
# Verifies required command-line tools and Python packages are installed.
# Run this first — fail fast rather than halfway through a download.

set -euo pipefail

echo "=== Phase 1 setup check ==="
echo

MISSING=()

check_tool () {
    local tool="$1"
    local install_hint="$2"
    if command -v "$tool" >/dev/null 2>&1; then
        printf "  [OK]   %s  (%s)\n" "$tool" "$(command -v "$tool")"
    else
        printf "  [MISS] %s  — install: %s\n" "$tool" "$install_hint"
        MISSING+=("$tool")
    fi
}

echo "Command-line tools:"
check_tool curl       "apt install curl"
check_tool gunzip     "apt install gzip"
check_tool bedtools   "conda install -c bioconda bedtools   OR   apt install bedtools"
check_tool liftOver   "conda install -c bioconda ucsc-liftover   OR   download from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
check_tool bigWigAverageOverBed "conda install -c bioconda ucsc-bigwigaverageoverbed   (optional but used in Phase 2)"

echo
echo "Python packages:"
python3 - <<'PY'
import importlib, sys
pkgs = ["pandas", "numpy"]
missing = []
for p in pkgs:
    try:
        m = importlib.import_module(p)
        print(f"  [OK]   {p:10s} {getattr(m, '__version__', '?')}")
    except ImportError:
        print(f"  [MISS] {p:10s} — install: pip install {p}")
        missing.append(p)
sys.exit(1 if missing else 0)
PY

if [ ${#MISSING[@]} -gt 0 ]; then
    echo
    echo "Missing tools: ${MISSING[*]}"
    echo "Install them and re-run this check."
    exit 1
fi

echo
echo "All required tools found. Ready to proceed."
