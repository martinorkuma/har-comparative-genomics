"""
_build_matched_cnes.py

Sample a control set of conserved non-accelerated noncoding elements (CNEs) that is
matched to the HARs on length and mean phastCons score.

Strategy:
    1. Start from UCSC phastCons 100-way pre-called elements (the candidate pool).
    2. Drop any candidate that:
        - overlaps a HAR (with a buffer, e.g. ±10 kb), or
        - overlaps a coding exon (we want noncoding controls).
    3. For each HAR, compute (length, mean_phastCons).
    4. Sample CNEs that are within length and conservation tolerances of the HAR set,
       at a `match_ratio` (default 3:1).

The output BED has the same columns as `hars.hg38.bed`, with names like CNE_000123.

Run via `bash scripts/01_acquire.sh` — not directly. Reads `config.yaml`.
"""
from __future__ import annotations

import sys
import yaml
import random
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools
from pathlib import Path

CFG = yaml.safe_load(open("config.yaml"))


def mean_score(bw, chrom: str, start: int, end: int) -> float:
    try:
        v = bw.stats(chrom, start, end, type="mean", exact=True)[0]
        return float(v) if v is not None else np.nan
    except RuntimeError:
        return np.nan


def main() -> None:
    rng = random.Random(CFG["cnes"]["random_seed"])
    np.random.seed(CFG["cnes"]["random_seed"])

    hars_path     = Path("data/processed/hars.hg38.bed")
    elements_path = Path("data/raw/hg38.phastConsElements100way.bed")
    bw_path       = Path("data/raw/hg38.phastCons100way.bw")
    gtf_gz        = Path("data/raw/gencode.v45.basic.annotation.gtf.gz")
    out_path      = Path("data/processed/cnes.hg38.bed")

    # --- coding exon BED (to exclude exonic phastCons elements) -----------
    exon_bed = Path("data/processed/_coding_exons.bed")
    if not exon_bed.exists():
        import gzip, re
        print("[cnes] extracting coding exons from GENCODE...")
        with gzip.open(gtf_gz, "rt") as fh, exon_bed.open("w") as out:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[2] != "CDS":
                    continue
                out.write(f"{f[0]}\t{int(f[3])-1}\t{f[4]}\n")
        # sort + merge
        pybedtools.BedTool(str(exon_bed)).sort().merge().saveas(str(exon_bed))

    hars       = pybedtools.BedTool(str(hars_path))
    elements   = pybedtools.BedTool(str(elements_path))
    exons      = pybedtools.BedTool(str(exon_bed))
    buffer_bp  = CFG["cnes"]["buffer_bp"]

    # --- candidate pool: phastCons elements minus HARs(±buffer) minus coding exons
    print(f"[cnes] candidate pool starts with {elements.count()} phastCons elements")
    hars_buffered = hars.slop(b=buffer_bp, genome="hg38")
    candidates = (
        elements
        .intersect(hars_buffered, v=True)   # not overlapping HAR ± buffer
        .intersect(exons,         v=True)   # noncoding only
    )
    print(f"[cnes] after dropping HAR-proximal and exonic: {candidates.count()} candidates")

    # --- score every HAR and every candidate with mean phastCons ----------
    bw = pyBigWig.open(str(bw_path))

    def score_bed(bt: pybedtools.BedTool) -> pd.DataFrame:
        rows = []
        for iv in bt:
            chrom, start, end = iv.chrom, iv.start, iv.end
            length = end - start
            mean   = mean_score(bw, chrom, start, end)
            rows.append((chrom, start, end, iv.name, length, mean))
        return pd.DataFrame(rows, columns=["chrom","start","end","name","length","mean_phastcons"])

    print("[cnes] scoring HARs ...")
    har_df = score_bed(hars).dropna(subset=["mean_phastcons"])
    print("[cnes] scoring candidate CNEs ...")
    cand_df = score_bed(candidates).dropna(subset=["mean_phastcons"])

    # --- match each HAR to k candidates within tolerance ------------------
    n_target_per_har = CFG["cnes"]["match_ratio"]
    len_tol  = CFG["cnes"]["length_tolerance_pct"] / 100.0
    cons_tol = CFG["cnes"]["conservation_tolerance"]

    # Sort candidates by length for fast lookup
    cand_df = cand_df.sort_values("length").reset_index(drop=True)
    used = set()

    picks = []
    for _, h in har_df.iterrows():
        lmin = h["length"] * (1 - len_tol)
        lmax = h["length"] * (1 + len_tol)
        cmin = h["mean_phastcons"] - cons_tol
        cmax = h["mean_phastcons"] + cons_tol
        # boolean mask of eligible candidates
        mask = (
            (cand_df["length"] >= lmin) & (cand_df["length"] <= lmax) &
            (cand_df["mean_phastcons"] >= cmin) & (cand_df["mean_phastcons"] <= cmax) &
            (~cand_df.index.isin(used))
        )
        eligible = cand_df.index[mask].tolist()
        rng.shuffle(eligible)
        for idx in eligible[:n_target_per_har]:
            used.add(idx)
            picks.append(idx)

    sampled = cand_df.loc[picks].copy()
    sampled["name"] = [f"CNE_{i:06d}" for i in range(len(sampled))]
    print(f"[cnes] sampled {len(sampled)} CNEs ({len(sampled)/len(har_df):.2f} per HAR)")

    # write BED6
    out = sampled[["chrom","start","end","name"]].copy()
    out["score"]  = 0
    out["strand"] = "+"
    out = out.sort_values(["chrom","start"])
    out.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"[cnes] wrote -> {out_path}")

    # write a quick QC summary
    qc = pd.DataFrame({
        "set":     ["HAR", "CNE"],
        "n":       [len(har_df), len(sampled)],
        "len_med": [har_df["length"].median(), sampled["length"].median()],
        "len_iqr": [har_df["length"].quantile(0.75)-har_df["length"].quantile(0.25),
                    sampled["length"].quantile(0.75)-sampled["length"].quantile(0.25)],
        "phastcons_mean": [har_df["mean_phastcons"].mean(), sampled["mean_phastcons"].mean()],
    })
    qc.to_csv("outputs/tables/cne_matching_qc.tsv", sep="\t", index=False)
    print("[cnes] QC -> outputs/tables/cne_matching_qc.tsv")


if __name__ == "__main__":
    sys.exit(main())
