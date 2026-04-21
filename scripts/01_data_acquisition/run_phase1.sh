#!/usr/bin/env bash
# run_phase1.sh
# Executes Phase 1 of the HAR vs CNE comparative genomics pipeline:
#   step 0 — verify tooling
#   step 1 — download HAR BED files (Pollard lab)
#   step 2 — download phastCons, GENCODE, liftOver chain
#   step 3 — liftOver HARs hg18 -> hg38
#   step 4 — build length-matched CNE controls
#
# Run from anywhere:   bash scripts/01_data_acquisition/run_phase1.sh
# Or from the phase dir: ./run_phase1.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "================================================="
echo " HAR project — Phase 1 (data acquisition)"
echo " Project root: $PROJECT_ROOT"
echo "================================================="

bash "$SCRIPT_DIR/00_setup_check.sh"
bash "$SCRIPT_DIR/01_download_hars.sh"
bash "$SCRIPT_DIR/02_download_references.sh"
bash "$SCRIPT_DIR/03_liftover_hars.sh"

python3 "$SCRIPT_DIR/04_build_matched_cnes.py" \
    --hars       "$PROJECT_ROOT/data/processed/hars.hg38.bed" \
    --phastcons  "$PROJECT_ROOT/data/raw/phastConsElements100way.hg38.bed" \
    --exons      "$PROJECT_ROOT/data/raw/gencode.v45.coding_exons.hg38.bed" \
    --out        "$PROJECT_ROOT/data/processed/cnes.hg38.bed" \
    --ratio 1 --length-tol 0.2 --seed 42 \
    2>&1 | tee "$PROJECT_ROOT/logs/04_build_matched_cnes.log"

echo
echo "================================================="
echo " Phase 1 complete. Outputs:"
echo "   $PROJECT_ROOT/data/processed/hars.hg38.bed"
echo "   $PROJECT_ROOT/data/processed/cnes.hg38.bed"
echo "================================================="
