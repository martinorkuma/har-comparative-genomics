"""
05_sensitivity_brain_tss.py

Sensitivity analysis for the reversed brain-TSS proximity result.

The headline finding from 04_interpret.py is that HARs sit FARTHER from
brain-expressed TSSes than the matched CNE controls (Spearman +0.578 on
dist_to_nearest_brain_tss). This script tests whether that reversal is:
    (a) a property of HARs, or
    (b) an artifact of matching CNEs on conservation, since conserved
        noncoding elements concentrate near genes.

Generates a third comparison set: random length-matched genomic intervals,
matched on length only (not conservation) and excluded from HARs and
coding exons so they remain noncoding. Compares distance-to-nearest-brain-TSS
distributions across HAR, matched CNE, and random sets.

Pre-decided interpretations (committed BEFORE looking at the result):
    Random > Matched CNE > HAR     -> matching artifact (strongest story)
    Random ~ Matched CNE > HAR     -> reversal is real, not an artifact
    Random ~ HAR < Matched CNE     -> matched CNEs are anomalously proximal

Outputs:
    outputs/figures/sensitivity_brain_tss.png
    outputs/tables/sensitivity_brain_tss.tsv

Run AFTER the main pipeline has finished:
    conda activate har-ml
    python scripts/05_sensitivity_brain_tss.py
"""
from __future__ import annotations

import sys
import yaml
import numpy as np
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu

CFG = yaml.safe_load(open("config.yaml"))
SEED = CFG["modeling"]["random_seed"]


def distances_to_nearest(query_bed: pybedtools.BedTool,
                         target_bed: pybedtools.BedTool) -> list[float]:
    """Absolute distance (bp) from each query to its nearest target. 0 = overlap."""
    closest = query_bed.sort().closest(target_bed.sort(), d=True, t="first")
    out = []
    for iv in closest:
        try:
            d = int(iv[-1])
            out.append(abs(d) if d >= 0 else np.nan)
        except (ValueError, IndexError):
            out.append(np.nan)
    return out


def build_chrom_sizes(fai_path: str, out_path: str) -> None:
    """Write a chrom.sizes file restricted to canonical chromosomes (chr1..22,X,Y)."""
    canonical = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    with open(fai_path) as fai, open(out_path, "w") as out:
        for line in fai:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in canonical:
                out.write(f"{parts[0]}\t{parts[1]}\n")


def main() -> None:
    np.random.seed(SEED)

    hars_path      = "data/processed/hars.hg38.bed"
    cnes_path      = "data/processed/cnes.hg38.bed"
    brain_tss_path = "data/processed/brain_tss.hg38.bed"
    exon_bed_path  = "data/processed/_coding_exons.bed"
    fai_path       = "data/raw/hg38.fa.fai"

    for p in (hars_path, cnes_path, brain_tss_path, exon_bed_path, fai_path):
        if not Path(p).exists():
            print(f"[sensitivity] ERROR: missing input {p}", file=sys.stderr)
            print("[sensitivity] run the main pipeline first.", file=sys.stderr)
            sys.exit(1)

    # --- chromosome sizes (canonical only) --------------------------------
    chrom_sizes_path = "data/processed/_hg38.canonical.chrom.sizes"
    if not Path(chrom_sizes_path).exists():
        build_chrom_sizes(fai_path, chrom_sizes_path)

    # --- inputs -----------------------------------------------------------
    hars      = pybedtools.BedTool(hars_path)
    cnes      = pybedtools.BedTool(cnes_path)
    brain_tss = pybedtools.BedTool(brain_tss_path)

    # --- exclusion bed: HARs + coding exons -------------------------------
    # keeps random intervals from overlapping the positive set or coding sequence
    print("[sensitivity] building exclusion regions (HARs + coding exons) ...")
    excl_path = "data/processed/_random_excl.bed"
    pybedtools.BedTool(hars_path) \
        .cat(exon_bed_path, postmerge=True) \
        .sort() \
        .saveas(excl_path)

    # --- random length-matched intervals ----------------------------------
    print(f"[sensitivity] shuffling {len(hars)} length-matched random intervals ...")
    random_bed = hars.shuffle(
        g=chrom_sizes_path,
        excl=excl_path,
        chrom=False,          # any chromosome, weighted by chromosome size
        noOverlapping=True,   # random intervals don't overlap each other
        seed=SEED,
    )
    n_har, n_cne, n_rand = len(hars), len(cnes), random_bed.count()
    print(f"[sensitivity]   HAR={n_har}  matched_CNE={n_cne}  random={n_rand}")

    # --- distances to brain TSS -------------------------------------------
    print("[sensitivity] computing distance to nearest brain-expressed TSS ...")
    har_d  = distances_to_nearest(hars,       brain_tss)
    cne_d  = distances_to_nearest(cnes,       brain_tss)
    rand_d = distances_to_nearest(random_bed, brain_tss)

    df = pd.DataFrame({
        "set": (["HAR"] * len(har_d)
                + ["Matched CNE"] * len(cne_d)
                + ["Random length-matched"] * len(rand_d)),
        "dist_bp": har_d + cne_d + rand_d,
    }).dropna(subset=["dist_bp"])
    df["log_dist"] = np.log1p(df["dist_bp"])

    # --- summary stats ----------------------------------------------------
    summ = (df.groupby("set")["dist_bp"]
              .describe(percentiles=[0.25, 0.5, 0.75])
              [["count", "mean", "25%", "50%", "75%"]]
              .rename(columns={"50%": "median", "25%": "q25", "75%": "q75"})
              .round(0)
              .astype({"count": int}))
    print("\nDistance to nearest brain TSS (bp):")
    print(summ.to_string(), "\n")

    # --- Mann-Whitney U pairwise ------------------------------------------
    har_arr  = np.asarray(har_d,  dtype=float); har_arr  = har_arr[~np.isnan(har_arr)]
    cne_arr  = np.asarray(cne_d,  dtype=float); cne_arr  = cne_arr[~np.isnan(cne_arr)]
    rand_arr = np.asarray(rand_d, dtype=float); rand_arr = rand_arr[~np.isnan(rand_arr)]

    u1, p1 = mannwhitneyu(har_arr, cne_arr,  alternative="two-sided")
    u2, p2 = mannwhitneyu(har_arr, rand_arr, alternative="two-sided")
    u3, p3 = mannwhitneyu(cne_arr, rand_arr, alternative="two-sided")
    print("Mann-Whitney U (two-sided):")
    print(f"  HAR vs Matched CNE:    U={u1:.4g}  p={p1:.4g}")
    print(f"  HAR vs Random:         U={u2:.4g}  p={p2:.4g}")
    print(f"  Matched CNE vs Random: U={u3:.4g}  p={p3:.4g}\n")

    # --- pre-decided interpretations (printed for the runner) -------------
    print("Pre-decided interpretations (compare medians above):")
    print("  Random > Matched CNE > HAR  -> matching artifact (strongest story)")
    print("  Random ~ Matched CNE > HAR  -> reversal is real, not an artifact")
    print("  Random ~ HAR < Matched CNE  -> matched CNEs anomalously proximal")
    print()

    # --- write the table ---------------------------------------------------
    summ.to_csv("outputs/tables/sensitivity_brain_tss.tsv", sep="\t")
    with open("outputs/tables/sensitivity_brain_tss.tsv", "a") as fh:
        fh.write("\n# Mann-Whitney U (two-sided), distance to nearest brain TSS:\n")
        fh.write(f"# HAR vs Matched CNE:    U={u1:.4g}  p={p1:.4g}\n")
        fh.write(f"# HAR vs Random:         U={u2:.4g}  p={p2:.4g}\n")
        fh.write(f"# Matched CNE vs Random: U={u3:.4g}  p={p3:.4g}\n")
    print("[sensitivity] -> outputs/tables/sensitivity_brain_tss.tsv")

    # --- plot --------------------------------------------------------------
    print("[sensitivity] plotting three-way comparison ...")
    plt.rcParams.update({"font.size": 11,
                         "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    order = ["HAR", "Matched CNE", "Random length-matched"]
    palette = {
        "HAR": "#cc4d4d",
        "Matched CNE": "#888888",
        "Random length-matched": "#4d7fcc",
    }
    sns.violinplot(data=df, x="set", y="log_dist", order=order,
                   palette=palette, inner="quartile", cut=0, ax=ax)
    # overlay medians as solid black bars
    medians = df.groupby("set")["log_dist"].median().reindex(order)
    for i, m in enumerate(medians.values):
        ax.plot([i - 0.18, i + 0.18], [m, m], color="black", lw=2)
    ax.set_xlabel("")
    ax.set_ylabel("log(1 + distance to nearest brain-expressed TSS, bp)")
    ax.set_title("Sensitivity: brain-TSS proximity by element set", pad=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig("outputs/figures/sensitivity_brain_tss.png",
                dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[sensitivity] -> outputs/figures/sensitivity_brain_tss.png")


if __name__ == "__main__":
    sys.exit(main())
