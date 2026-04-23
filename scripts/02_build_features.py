"""
02_build_features.py

Compute the 7 features per element (HAR or CNE) and write a tidy table:
    data/processed/features.tsv

Columns:
    chrom, start, end, name, label
    gc_content
    length
    mean_phastcons
    dist_to_nearest_tss
    dist_to_nearest_brain_tss
    overlaps_ccre              (0/1)
    overlaps_fetal_brain_enh   (0/1)

label = 1 for HAR, 0 for CNE.
"""
from __future__ import annotations

import sys
import yaml
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools
from pathlib import Path
from pyfaidx import Fasta

CFG = yaml.safe_load(open("config.yaml"))


def gc_content(seq: str) -> float:
    seq = seq.upper()
    n = sum(1 for b in seq if b in "ACGT")
    if n == 0:
        return np.nan
    g = sum(1 for b in seq if b in "GC")
    return g / n


def mean_phastcons(bw, chrom: str, start: int, end: int) -> float:
    try:
        v = bw.stats(chrom, start, end, type="mean", exact=True)[0]
        return float(v) if v is not None else np.nan
    except RuntimeError:
        return np.nan


def load_bed_with_label(path: str, label: int) -> pd.DataFrame:
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom","start","end","name","score","strand"],
        usecols=[0,1,2,3],
        dtype={"chrom":str, "start":int, "end":int, "name":str},
    )
    df["label"] = label
    return df


def signed_distance_to_nearest(bed_query: pybedtools.BedTool,
                               bed_target: pybedtools.BedTool) -> pd.Series:
    """For each interval in query, distance to nearest interval in target.
    Returns absolute distance in bp (0 if overlapping)."""
    closest = bed_query.sort().closest(bed_target.sort(), d=True, t="first")
    rows = []
    for iv in closest:
        # The d field is the last column (an integer).
        # Index from the end to avoid shifting if BED widths differ.
        try:
            d = int(iv[-1])
        except (ValueError, IndexError):
            d = -1
        rows.append((iv.name, abs(d) if d >= 0 else np.nan))
    return pd.DataFrame(rows, columns=["name","dist"]).set_index("name")["dist"]


def overlaps_any(bed_query: pybedtools.BedTool,
                 bed_target: pybedtools.BedTool) -> pd.Series:
    """Returns 1 if a query interval overlaps any target interval, else 0."""
    overlapping_names = {iv.name for iv in bed_query.intersect(bed_target, u=True)}
    all_names = [iv.name for iv in bed_query]
    return pd.Series(
        [1 if n in overlapping_names else 0 for n in all_names],
        index=all_names, name="overlaps",
    )


def main() -> None:
    print("[features] loading HARs and CNEs ...")
    hars = load_bed_with_label("data/processed/hars.hg38.bed", label=1)
    cnes = load_bed_with_label("data/processed/cnes.hg38.bed", label=0)
    elems = pd.concat([hars, cnes], ignore_index=True)
    print(f"[features]   {len(hars)} HARs + {len(cnes)} CNEs = {len(elems)} elements")

    elems["length"] = elems["end"] - elems["start"]

    # --- GC content -------------------------------------------------------
    print("[features] computing GC content ...")
    fa = Fasta("data/raw/hg38.fa", as_raw=True, sequence_always_upper=True)
    elems["gc_content"] = [
        gc_content(fa[r["chrom"]][r["start"]:r["end"]])
        for _, r in elems.iterrows()
    ]

    # --- mean phastCons ---------------------------------------------------
    print("[features] computing mean phastCons ...")
    bw = pyBigWig.open("data/raw/hg38.phastCons100way.bw")
    elems["mean_phastcons"] = [
        mean_phastcons(bw, r["chrom"], r["start"], r["end"])
        for _, r in elems.iterrows()
    ]

    # --- Distances and overlaps via BedTools ------------------------------
    # write a temp BED for the combined element set, named by `name`
    tmp_bed = "data/processed/_features_query.bed"
    elems[["chrom","start","end","name"]].assign(score=0, strand="+") \
        .to_csv(tmp_bed, sep="\t", header=False, index=False)
    query = pybedtools.BedTool(tmp_bed).sort()

    print("[features] distance to nearest TSS (any GENCODE gene) ...")
    all_tss   = pybedtools.BedTool("data/processed/all_tss.hg38.bed").sort()
    elems = elems.merge(
        signed_distance_to_nearest(query, all_tss).rename("dist_to_nearest_tss"),
        left_on="name", right_index=True,
    )

    print("[features] distance to nearest brain-expressed TSS ...")
    brain_tss = pybedtools.BedTool("data/processed/brain_tss.hg38.bed").sort()
    elems = elems.merge(
        signed_distance_to_nearest(query, brain_tss).rename("dist_to_nearest_brain_tss"),
        left_on="name", right_index=True,
    )

    print("[features] overlap with ENCODE cCREs ...")
    ccre = pybedtools.BedTool("data/processed/encode_ccre.hg38.bed").sort()
    elems = elems.merge(
        overlaps_any(query, ccre).rename("overlaps_ccre").to_frame(),
        left_on="name", right_index=True,
    )

    print("[features] overlap with fetal-brain active enhancers ...")
    fbe = pybedtools.BedTool("data/processed/fetal_brain_enhancers.hg38.bed").sort()
    elems = elems.merge(
        overlaps_any(query, fbe).rename("overlaps_fetal_brain_enh").to_frame(),
        left_on="name", right_index=True,
    )

    # --- finalize ---------------------------------------------------------
    feature_cols = [
        "gc_content", "length", "mean_phastcons",
        "dist_to_nearest_tss", "dist_to_nearest_brain_tss",
        "overlaps_ccre", "overlaps_fetal_brain_enh",
    ]
    out = elems[["chrom","start","end","name","label", *feature_cols]]
    out = out.dropna(subset=feature_cols)  # drop rows where any feature is NaN

    n_har_kept = (out["label"]==1).sum()
    n_cne_kept = (out["label"]==0).sum()
    n_har_drop = len(hars) - n_har_kept
    n_cne_drop = len(cnes) - n_cne_kept
    print(f"[features] kept {n_har_kept} HARs, {n_cne_kept} CNEs "
          f"(dropped {n_har_drop}+{n_cne_drop} due to missing features)")

    out.to_csv("data/processed/features.tsv", sep="\t", index=False)
    print(f"[features] wrote -> data/processed/features.tsv")

    # quick summary table for the poster
    summ = (out.groupby("label")[feature_cols]
              .agg(["median","mean"])
              .round(4))
    summ.to_csv("outputs/tables/feature_summary.tsv", sep="\t")
    print("[features] summary -> outputs/tables/feature_summary.tsv")


if __name__ == "__main__":
    sys.exit(main())
