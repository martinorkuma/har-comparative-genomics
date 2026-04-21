#!/usr/bin/env bash
# 03_liftover_hars.sh
# Lifts HAR coordinates from hg18 -> hg38 using UCSC liftOver.
#
# The Pollard-lab HAR BEDs are in hg18. Modern annotation resources
# (phastCons100way, ENCODE cCREs, GENCODE v45) are in hg38. We must align
# coordinate systems before any overlap analysis.
#
# Outputs:
#   data/processed/hars.primate.hg38.bed    (successfully lifted primate HARs)
#   data/processed/hars.primate.hg38.unmapped  (regions that failed to map)
#   data/processed/hars.mammal.hg38.bed     (mammal set)
#   data/processed/hars.mammal.hg38.unmapped
#   data/processed/hars.hg38.bed            (combined, de-duplicated union)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DIR="$PROJECT_ROOT/data/raw"
PROC_DIR="$PROJECT_ROOT/data/processed"
LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$PROC_DIR" "$LOG_DIR"

LOG="$LOG_DIR/03_liftover_hars.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) — 03_liftover_hars.sh ==="

CHAIN_GZ="$RAW_DIR/hg18ToHg38.over.chain.gz"
CHAIN="$RAW_DIR/hg18ToHg38.over.chain"

if [ ! -s "$CHAIN" ]; then
    if [ ! -s "$CHAIN_GZ" ]; then
        echo "ERROR: chain file missing. Run 02_download_references.sh first."
        exit 1
    fi
    echo "  [cvt]  Decompressing liftOver chain ..."
    gunzip -kc "$CHAIN_GZ" > "$CHAIN"
fi

liftover_one () {
    local label="$1"       # "primate" or "mammal"
    local hg18_bed="$2"    # input BED (hg18)
    if [ ! -s "$hg18_bed" ]; then
        echo "  [skip] $hg18_bed missing — skipping $label set"
        return 0
    fi
    local out_bed="$PROC_DIR/hars.${label}.hg38.bed"
    local out_unmapped="$PROC_DIR/hars.${label}.hg38.unmapped"
    local tagged="$PROC_DIR/.tmp.${label}.tagged.bed"

    # Assign stable IDs to input so we can trace lifted regions back.
    # Pollard BEDs may have fewer than 4 columns; we force BED4 with a label.
    awk -v pfx="HAR_${label}_" 'BEGIN{OFS="\t"} !/^track|^browser|^#/ && NF>=3 {
        id = pfx sprintf("%05d", ++i)
        print $1, $2, $3, id
    }' "$hg18_bed" > "$tagged"

    echo "  [lft]  $label: $(wc -l < "$tagged") input HARs"
    liftOver -minMatch=0.5 "$tagged" "$CHAIN" "$out_bed" "$out_unmapped"
    local mapped=$(wc -l < "$out_bed")
    local unmapped=$(grep -vc '^#' "$out_unmapped" 2>/dev/null || echo 0)
    echo "  [done] $label: $mapped mapped, $unmapped unmapped -> $out_bed"
    rm -f "$tagged"
}

liftover_one primate "$RAW_DIR/2xPrimateHARs.hg18.bed"
liftover_one mammal  "$RAW_DIR/2xMammalHARs.hg18.bed"

# ---- Build combined union HAR set ----
# Merge primate + mammal sets. If they overlap, bedtools merge will join them;
# we keep the primate ID (arbitrary but deterministic — primate is larger).
COMBINED="$PROC_DIR/hars.hg38.bed"
{
    [ -s "$PROC_DIR/hars.primate.hg38.bed" ] && cat "$PROC_DIR/hars.primate.hg38.bed"
    [ -s "$PROC_DIR/hars.mammal.hg38.bed"  ] && cat "$PROC_DIR/hars.mammal.hg38.bed"
} | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
  | bedtools merge -i - -c 4 -o first > "$COMBINED"

echo "  [done] combined union -> $COMBINED ($(wc -l < "$COMBINED") unique HARs)"
echo "=== Done ==="
