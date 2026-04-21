#!/usr/bin/env bash
# 02_download_references.sh
# Downloads reference annotation data needed for CNE control construction:
#
#   1. phastConsElements100way (hg38): conserved elements from 100-way
#      vertebrate alignment — the pool we'll draw CNE controls from.
#   2. GENCODE basic annotation (hg38): to identify coding exons, which we
#      exclude when building non-coding CNEs.
#   3. hg18ToHg38 liftOver chain: needed to bring HARs into hg38.
#
# All files are stored under data/raw/. The script is idempotent — re-running
# will skip files that already exist and are non-empty.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RAW_DIR="$PROJECT_ROOT/data/raw"
LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$RAW_DIR" "$LOG_DIR"

LOG="$LOG_DIR/02_download_references.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) — 02_download_references.sh ==="

# ---- URLs (UCSC goldenPath + GENCODE FTP) ----
PHASTCONS_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/phastConsElements100way.txt.gz"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz"
# GENCODE v45 basic annotation (stable, widely used; adjust version if needed)
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz"

PHASTCONS_GZ="$RAW_DIR/phastConsElements100way.txt.gz"
CHAIN_GZ="$RAW_DIR/hg18ToHg38.over.chain.gz"
GENCODE_GZ="$RAW_DIR/gencode.v45.basic.annotation.gtf.gz"

download () {
    local url="$1" out="$2"
    if [ -s "$out" ]; then
        echo "  [skip] $out ($(du -h "$out" | cut -f1))"
        return 0
    fi
    echo "  [get]  $url"
    curl -fLsS --retry 3 --retry-delay 2 -o "$out" "$url"
    echo "  [done] $out ($(du -h "$out" | cut -f1))"
}

download "$PHASTCONS_URL" "$PHASTCONS_GZ"
download "$CHAIN_URL"     "$CHAIN_GZ"
download "$GENCODE_URL"   "$GENCODE_GZ"

# ---- Convert phastCons UCSC MySQL dump to BED ----
# Columns in phastConsElements100way.txt.gz:
#   bin  chrom  chromStart  chromEnd  name  score  lod
# We output BED6: chrom, start, end, name, lod(score), strand(.)
# `lod` is the log-odds conservation score — preserved in column 5 so we can
# length- AND conservation-match CNEs to HARs downstream.
PHASTCONS_BED="$RAW_DIR/phastConsElements100way.hg38.bed"
if [ ! -s "$PHASTCONS_BED" ]; then
    echo "  [cvt]  phastCons text -> BED6 ..."
    zcat "$PHASTCONS_GZ" | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $7, "."}' \
        | sort -k1,1 -k2,2n > "$PHASTCONS_BED"
    echo "  [done] $PHASTCONS_BED ($(wc -l < "$PHASTCONS_BED") elements)"
else
    echo "  [skip] $PHASTCONS_BED ($(wc -l < "$PHASTCONS_BED") elements)"
fi

# ---- Extract coding exons from GENCODE GTF -> BED ----
# Keep only features with feature_type == "exon" AND transcript_type == "protein_coding".
# We merge overlapping exons with bedtools downstream in CNE construction.
EXONS_BED="$RAW_DIR/gencode.v45.coding_exons.hg38.bed"
if [ ! -s "$EXONS_BED" ]; then
    echo "  [cvt]  GENCODE GTF -> coding exons BED ..."
    zcat "$GENCODE_GZ" \
        | awk 'BEGIN{OFS="\t"} !/^#/ && $3=="exon" && /transcript_type "protein_coding"/ {
                 # GTF is 1-based, inclusive; BED is 0-based, half-open
                 print $1, $4-1, $5, ".", ".", $7
               }' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - > "$EXONS_BED"
    echo "  [done] $EXONS_BED ($(wc -l < "$EXONS_BED") merged exon intervals)"
else
    echo "  [skip] $EXONS_BED ($(wc -l < "$EXONS_BED") merged exon intervals)"
fi

echo "=== Done: references in $RAW_DIR ==="
