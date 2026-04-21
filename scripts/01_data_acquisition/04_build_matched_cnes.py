#!/usr/bin/env python3
"""
04_build_matched_cnes.py

Construct a matched control set of Conserved Non-accelerated noncoding
Elements (CNEs) for the HAR positive set.

Pipeline:
  1. Load phastCons100way elements (hg38) — the pool of candidate conserved
     elements.
  2. Filter out elements that overlap protein-coding exons (we want NON-
     coding conserved elements, matching HARs which are predominantly
     noncoding).
  3. Filter out elements that overlap any HAR (controls must not be HARs).
  4. Restrict both sets to standard autosomes + chrX (drop Y, MT, alts,
     unplaced contigs).
  5. For each HAR, sample N control CNEs whose length falls within a
     tolerance band of the HAR's length. Sampling is without replacement
     across the whole run so CNEs are never reused.

The output BED4 carries a name field like "CNE_for_HAR_primate_00042" so
each control is traceable to the HAR it was matched against. This is
useful later if we do paired analyses or stratified interpretation.

Usage:
  python 04_build_matched_cnes.py \
      --hars       data/processed/hars.hg38.bed \
      --phastcons  data/raw/phastConsElements100way.hg38.bed \
      --exons      data/raw/gencode.v45.coding_exons.hg38.bed \
      --out        data/processed/cnes.hg38.bed \
      --ratio 1 --length-tol 0.2 --seed 42
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX"}

BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )


def run_bedtools_subtract(a: Path, b: Path) -> Path:
    """Return a temp file with (a - b). Keeps full a-records; does NOT split
    them at overlap boundaries (-A flag)."""
    tmp = Path(tempfile.mkstemp(suffix=".bed")[1])
    cmd = ["bedtools", "subtract", "-A", "-a", str(a), "-b", str(b)]
    with tmp.open("w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)
    return tmp


def load_bed(path: Path, min_cols: int = 3) -> pd.DataFrame:
    """Load a BED file, filling missing optional columns with '.'."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     dtype={0: str}, low_memory=False)
    ncols = df.shape[1]
    if ncols < min_cols:
        raise ValueError(f"{path} has {ncols} columns, need at least {min_cols}")
    # Assign standard BED col names for the cols we have
    df.columns = BED_COLS[:ncols] + list(df.columns[len(BED_COLS):])
    # Ensure numeric coords
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def filter_standard_chroms(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)
    df = df[df["chrom"].isin(STANDARD_CHROMS)].copy()
    logging.info("  chrom filter: %d -> %d (dropped %d on non-standard contigs)",
                 before, len(df), before - len(df))
    return df


def match_cnes_to_hars(hars: pd.DataFrame,
                       candidates: pd.DataFrame,
                       ratio: int,
                       length_tol: float,
                       seed: int) -> pd.DataFrame:
    """
    For each HAR, draw `ratio` CNEs from candidates with length within
    length_tol fractional deviation (e.g. 0.2 = ±20 %). Sample without
    replacement across the whole run.
    """
    rng = np.random.default_rng(seed)
    hars = hars.copy()
    hars["length"] = hars["end"] - hars["start"]
    candidates = candidates.copy()
    candidates["length"] = candidates["end"] - candidates["start"]

    # Sort candidates by length once for fast windowed lookup via searchsorted
    candidates = candidates.sort_values("length").reset_index(drop=True)
    cand_lengths = candidates["length"].to_numpy()

    used = np.zeros(len(candidates), dtype=bool)
    picks = []
    skipped_no_candidates = 0

    # Iterate HARs in random order so early HARs don't monopolise the pool
    har_order = rng.permutation(len(hars))
    for hi in har_order:
        har = hars.iloc[hi]
        lo = har["length"] * (1 - length_tol)
        hi_ = har["length"] * (1 + length_tol)
        lo_idx = np.searchsorted(cand_lengths, lo, side="left")
        hi_idx = np.searchsorted(cand_lengths, hi_, side="right")
        pool_idx = np.arange(lo_idx, hi_idx)
        # Restrict to unused candidates
        pool_idx = pool_idx[~used[pool_idx]]
        if len(pool_idx) == 0:
            skipped_no_candidates += 1
            continue
        k = min(ratio, len(pool_idx))
        chosen = rng.choice(pool_idx, size=k, replace=False)
        used[chosen] = True
        for ci in chosen:
            cand = candidates.iloc[ci]
            picks.append({
                "chrom": cand["chrom"],
                "start": int(cand["start"]),
                "end":   int(cand["end"]),
                "name":  f"CNE_for_{har['name']}",
                "score": cand.get("score", "."),
                "strand": ".",
            })

    logging.info("  matched %d CNEs to %d HARs (ratio %d, tol ±%.0f%%)",
                 len(picks), len(hars) - skipped_no_candidates, ratio, length_tol * 100)
    if skipped_no_candidates:
        logging.warning("  %d HARs had no length-matched CNE available",
                        skipped_no_candidates)

    return pd.DataFrame(picks, columns=["chrom", "start", "end", "name", "score", "strand"])


def summarise(df: pd.DataFrame, label: str) -> None:
    if df.empty:
        logging.info("  %s: EMPTY", label)
        return
    ln = (df["end"] - df["start"])
    logging.info("  %s: n=%d  length median=%d  mean=%.0f  min=%d  max=%d",
                 label, len(df), int(ln.median()), ln.mean(), ln.min(), ln.max())


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--hars", required=True, type=Path)
    p.add_argument("--phastcons", required=True, type=Path)
    p.add_argument("--exons", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--ratio", type=int, default=1,
                   help="Number of CNE controls per HAR (default 1)")
    p.add_argument("--length-tol", type=float, default=0.2,
                   help="Fractional length tolerance, e.g. 0.2 = ±20%% (default 0.2)")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    setup_logging()
    logging.info("=== 04_build_matched_cnes.py ===")

    # 1. Load HARs
    logging.info("Loading HARs from %s ...", args.hars)
    hars = load_bed(args.hars, min_cols=4)
    hars = filter_standard_chroms(hars)
    summarise(hars, "HARs (standard chroms)")
    if hars.empty:
        logging.error("No HARs after chrom filter. Aborting.")
        return 1

    # 2. Subtract exons and HARs from phastCons pool using bedtools (faster than pandas)
    logging.info("Filtering phastCons against exons and HARs (bedtools subtract -A) ...")
    no_exons = run_bedtools_subtract(args.phastcons, args.exons)
    no_exons_no_hars = run_bedtools_subtract(no_exons, args.hars)
    try:
        candidates = load_bed(no_exons_no_hars, min_cols=3)
    finally:
        no_exons.unlink(missing_ok=True)
        no_exons_no_hars.unlink(missing_ok=True)
    candidates = filter_standard_chroms(candidates)
    summarise(candidates, "CNE candidates (conserved, non-coding, non-HAR)")
    if candidates.empty:
        logging.error("No CNE candidates remain. Aborting.")
        return 1

    # 3. Length-matched sampling
    logging.info("Matching CNE controls to HARs ...")
    cnes = match_cnes_to_hars(
        hars, candidates,
        ratio=args.ratio,
        length_tol=args.length_tol,
        seed=args.seed,
    )
    summarise(cnes, f"Selected CNEs ({args.ratio}:1)")

    # 4. Write output sorted
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cnes = cnes.sort_values(["chrom", "start"]).reset_index(drop=True)
    cnes.to_csv(args.out, sep="\t", header=False, index=False)
    logging.info("Wrote %d CNEs to %s", len(cnes), args.out)

    # 5. Print length distribution comparison
    import io
    buf = io.StringIO()
    pd.DataFrame({
        "HAR_length": (hars["end"] - hars["start"]).describe(),
        "CNE_length": (cnes["end"] - cnes["start"]).describe(),
    }).to_string(buf)
    logging.info("Length distributions:\n%s", buf.getvalue())

    return 0


if __name__ == "__main__":
    sys.exit(main())
