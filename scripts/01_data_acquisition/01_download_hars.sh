#!/usr/bin/env bash
# 01_download_hars.sh
# Downloads HAR coordinate BED files from the Pollard lab website.
#
# Sources:
#   Primate HARs: 1,930 elements, FDR < 0.1 (hg18 coordinates)
#     https://docpollard.org/research/human-acceleration-in-primate-conserved-elements/
#   Mammal HARs:    563 elements, FDR < 0.1 (hg18 coordinates)
#     https://docpollard.org/research/human-acceleration-in-mammal-conserved-elements/
#
# Both sets are in hg18 and must be lifted over to hg38 (see 03_liftover_hars.sh).

set -euo pipefail

# Resolve project root regardless of where script is invoked from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DIR="$PROJECT_ROOT/data/raw"
LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$RAW_DIR" "$LOG_DIR"

LOG="$LOG_DIR/01_download_hars.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) — 01_download_hars.sh ==="

PRIMATE_URL="https://docpollard.org/wordpress/wp-content/research/2xPrimateHARs.bed"
MAMMAL_URL="https://docpollard.org/wordpress/wp-content/research/2xMammalHARs.bed"

PRIMATE_OUT="$RAW_DIR/2xPrimateHARs.hg18.bed"
MAMMAL_OUT="$RAW_DIR/2xMammalHARs.hg18.bed"

download () {
    local url="$1"
    local out="$2"
    if [ -s "$out" ]; then
        echo "  [skip] $out already exists ($(wc -l < "$out") lines)"
        return 0
    fi
    echo "  [get]  $url"
    # -f: fail on HTTP errors; -L: follow redirects; --retry for flaky connections
    curl -fLsS --retry 3 --retry-delay 2 -o "$out" "$url"
    echo "  [done] $out ($(wc -l < "$out") lines)"
}

download "$PRIMATE_URL" "$PRIMATE_OUT"
download "$MAMMAL_URL"  "$MAMMAL_OUT"

# Sanity check: BED format and expected line counts (with small tolerance)
python3 - "$PRIMATE_OUT" "$MAMMAL_OUT" <<'PY'
import sys
files = sys.argv[1:]
expected = {"2xPrimateHARs.hg18.bed": 1930, "2xMammalHARs.hg18.bed": 563}
for f in files:
    name = f.rsplit("/", 1)[-1]
    with open(f) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith(("track", "#", "browser"))]
    exp = expected.get(name)
    ok_count = (exp is None) or (abs(len(lines) - exp) <= 5)
    # Check first line looks like BED (chrom, start, end)
    first = lines[0].rstrip("\n").split("\t")
    ok_format = len(first) >= 3 and first[0].startswith("chr") and first[1].isdigit() and first[2].isdigit()
    status = "OK" if (ok_count and ok_format) else "WARN"
    print(f"  [{status}] {name}: {len(lines)} lines, expected ~{exp}, first row: {first[:4]}")
PY

echo "=== Done: HAR BEDs in $RAW_DIR ==="
