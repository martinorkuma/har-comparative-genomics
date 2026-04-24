#!/usr/bin/env bash
# fix_line_endings.sh - strip CRLF from all scripts in the repo.
#
# The error you saw:
#     scripts/01_acquire.sh: line 17: ...: No such file or directory
# is the classic signature of Windows CRLF line endings in a bash script:
# bash reads '\r' as part of the command/heredoc delimiter and fails.
#
# Run this ONCE from the repo root, then replace 01_acquire.sh with the fixed
# version and run the pipeline again.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
cd "$ROOT"

echo "[fix] converting CRLF -> LF for shell and python scripts..."

# Prefer dos2unix if available, else sed.
if command -v dos2unix &>/dev/null; then
    find scripts -type f \( -name "*.sh" -o -name "*.py" \) -print0 \
        | xargs -0 dos2unix -q
    find . -maxdepth 1 -type f \( -name "*.sh" -o -name "*.yml" -o -name "*.yaml" \) -print0 \
        | xargs -0 dos2unix -q
else
    find scripts -type f \( -name "*.sh" -o -name "*.py" \) -exec sed -i 's/\r$//' {} +
    find . -maxdepth 1 -type f \( -name "*.sh" -o -name "*.yml" -o -name "*.yaml" \) \
        -exec sed -i 's/\r$//' {} +
fi

echo "[fix] done. verifying..."
for f in scripts/*.sh setup.sh; do
    if [[ -f "$f" ]] && grep -qU $'\r' "$f"; then
        echo "  STILL HAS CRLF: $f"
    else
        [[ -f "$f" ]] && echo "  ok: $f"
    fi
done

echo
echo "[fix] now run:"
echo "  bash -n scripts/01_acquire.sh   # syntax check"
echo "  bash scripts/run_all.sh         # re-run pipeline"
