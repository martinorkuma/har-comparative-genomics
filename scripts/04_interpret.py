"""
04_interpret.py

Generate the project's CENTRAL RESULT and the HARE5/FZD8 case study.

Outputs:
    outputs/figures/shap_summary.png          beeswarm SHAP summary (THE poster centerpiece)
    outputs/figures/shap_bar.png              mean |SHAP| ranking (cleaner version for talk)
    outputs/figures/feature_distributions.png HAR vs CNE per-feature distributions
    outputs/figures/hare5_case_study.png      HARE5's feature values vs HAR/CNE distributions
    outputs/tables/top_shap_features.tsv      ranked feature importances + direction of effect
    outputs/tables/hare5_features.tsv         HARE5's row from features.tsv
"""
from __future__ import annotations

import sys
import yaml
import joblib
import numpy as np
import pandas as pd
import shap
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

CFG = yaml.safe_load(open("config.yaml"))

FEATURE_COLS = [
    "gc_content", "length", "mean_phastcons",
    "dist_to_nearest_tss", "dist_to_nearest_brain_tss",
    "overlaps_ccre", "overlaps_fetal_brain_enh",
]
LOG_FEATURES = ["length", "dist_to_nearest_tss", "dist_to_nearest_brain_tss"]

# Pretty labels for plots
PRETTY = {
    "gc_content":               "GC content",
    "length":                   "log(length)",
    "mean_phastcons":           "Mean phastCons",
    "dist_to_nearest_tss":      "log(dist to nearest TSS)",
    "dist_to_nearest_brain_tss":"log(dist to nearest brain TSS)",
    "overlaps_ccre":            "Overlaps ENCODE cCRE",
    "overlaps_fetal_brain_enh": "Overlaps fetal-brain enhancer",
}


def transform(X: pd.DataFrame) -> pd.DataFrame:
    Xt = X.copy()
    for c in LOG_FEATURES:
        Xt[c] = np.log1p(Xt[c])
    return Xt


def liftover_hare5_interval() -> tuple[str, int, int] | None:
    """Lift Boyd 2015 HARE5 coordinates from hg19 to hg38."""
    import os
    import subprocess
    import tempfile

    chrom = CFG["case_study"]["hg19_chrom"]
    s_hg19 = CFG["case_study"]["hg19_start"]
    e_hg19 = CFG["case_study"]["hg19_end"]
    chain = "data/raw/hg19ToHg38.over.chain.gz"

    with tempfile.TemporaryDirectory() as td:
        in_bed = os.path.join(td, "hare5.hg19.bed")
        out_bed = os.path.join(td, "hare5.hg38.bed")
        unmapped = os.path.join(td, "hare5.unmapped.bed")
        with open(in_bed, "w") as fh:
            fh.write(f"{chrom}\t{s_hg19}\t{e_hg19}\tHARE5\n")
        subprocess.run(["liftOver", in_bed, chain, out_bed, unmapped], check=True)
        if os.path.getsize(out_bed) == 0:
            print("[interpret] WARNING: HARE5 failed liftOver; case study skipped.")
            return None
        with open(out_bed) as fh:
            f = fh.readline().rstrip("\n").split("\t")
            chrom38, s38, e38 = f[0], int(f[1]), int(f[2])

    print(
        "[interpret] HARE5/Boyd hg19 interval "
        f"{chrom}:{s_hg19:,}-{e_hg19:,} lifts to "
        f"{chrom38}:{s38:,}-{e38:,} (hg38)"
    )
    return chrom38, s38, e38


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def interval_gap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if overlap_bp(a_start, a_end, b_start, b_end) > 0:
        return 0
    if a_end <= b_start:
        return b_start - a_end
    return a_start - b_end


def build_external_hare5_row(chrom: str, start: int, end: int) -> pd.Series:
    """Compute HARE5 features when it is absent from the modeled HAR table."""
    import pyBigWig
    import pybedtools
    from pyfaidx import Fasta

    def gc(seq: str) -> float:
        seq = seq.upper()
        n = sum(1 for b in seq if b in "ACGT")
        return np.nan if n == 0 else sum(1 for b in seq if b in "GC") / n

    def bw_mean(path: str) -> float:
        with pyBigWig.open(path) as bw:
            v = bw.stats(chrom, start, end, type="mean", exact=True)[0]
        return float(v) if v is not None else np.nan

    def query_bed() -> pybedtools.BedTool:
        return pybedtools.BedTool(
            f"{chrom}\t{start}\t{end}\t{CFG['case_study']['name']}\n",
            from_string=True,
        ).sort()

    def nearest_distance(path: str) -> float:
        target = pybedtools.BedTool(path).sort()
        closest = query_bed().closest(target, d=True, t="first")
        try:
            return abs(int(next(iter(closest))[-1]))
        except (StopIteration, ValueError, IndexError):
            return np.nan

    def overlaps(path: str) -> int:
        target = pybedtools.BedTool(path).sort()
        return int(any(query_bed().intersect(target, u=True)))

    fa = Fasta("data/raw/hg38.fa", as_raw=True, sequence_always_upper=True)
    row = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "name": CFG["case_study"]["name"],
        "label": 1,
        "gc_content": gc(fa[chrom][start:end]),
        "length": end - start,
        "mean_phastcons": bw_mean("data/raw/hg38.phastCons100way.bw"),
        "dist_to_nearest_tss": nearest_distance("data/processed/all_tss.hg38.bed"),
        "dist_to_nearest_brain_tss": nearest_distance("data/processed/brain_tss.hg38.bed"),
        "overlaps_ccre": overlaps("data/processed/encode_ccre.hg38.bed"),
        "overlaps_fetal_brain_enh": overlaps("data/processed/fetal_brain_enhancers.hg38.bed"),
    }
    return pd.Series(row)


def find_hare5(df: pd.DataFrame) -> pd.Series | None:
    """Locate HARE5 in the feature table.

    Two attempts, in order:
      (a) liftOver HARE5's hg19 coordinates to hg38.
      (b) accept a name/coordinate match only if it overlaps the lifted interval.

    If HARE5 is absent from the modeled HAR table, compute an external HARE5 row
    from the same feature sources. This keeps the case-study figure honest while
    leaving the training/evaluation data unchanged.
    """
    lifted = liftover_hare5_interval()
    if lifted is None:
        return None
    chrom38, s38, e38 = lifted

    name_target = CFG["case_study"]["name"]
    hits = df[df["name"].astype(str).str.strip() == name_target]
    if len(hits):
        hit = hits.iloc[0]
        bp = overlap_bp(int(hit["start"]), int(hit["end"]), s38, e38)
        if hit["chrom"] == chrom38 and bp > 0:
            print(
                f"[interpret] case study: matched HAR by name '{name_target}' "
                f"at {hit['chrom']}:{int(hit['start']):,}-{int(hit['end']):,} "
                f"({bp:,} bp overlap with HARE5)"
            )
            return hit
        print(
            f"[interpret] WARNING: name '{name_target}' matched "
            f"{hit['chrom']}:{int(hit['start']):,}-{int(hit['end']):,}, "
            "but it does not overlap the Boyd HARE5 interval; ignoring name hit."
        )

    overlapping = df[
        (df["label"] == 1)
        & (df["chrom"] == chrom38)
        & (df["end"] > s38)
        & (df["start"] < e38)
    ].copy()
    if len(overlapping):
        overlapping["overlap_bp"] = overlapping.apply(
            lambda r: overlap_bp(int(r["start"]), int(r["end"]), s38, e38),
            axis=1,
        )
        hit = overlapping.sort_values("overlap_bp", ascending=False).iloc[0]
        print(
            f"[interpret] case study: matched HAR '{hit['name']}' by coordinate "
            f"at {hit['chrom']}:{int(hit['start']):,}-{int(hit['end']):,} "
            f"({int(hit['overlap_bp']):,} bp overlap with HARE5)"
        )
        return hit

    window = CFG["case_study"]["fallback_window_bp"]
    cand = df[(df["label"] == 1) & (df["chrom"] == chrom38)].copy()
    if cand.empty:
        print(f"[interpret] WARNING: no modeled HARs on {chrom38}.")
    else:
        cand["gap_bp"] = cand.apply(
            lambda r: interval_gap_bp(int(r["start"]), int(r["end"]), s38, e38),
            axis=1,
        )
        nearest = cand.sort_values("gap_bp").iloc[0]
        print(
            f"[interpret] HARE5 is absent from the modeled HAR table; nearest "
            f"HAR is '{nearest['name']}' at "
            f"{nearest['chrom']}:{int(nearest['start']):,}-{int(nearest['end']):,} "
            f"({int(nearest['gap_bp']):,} bp away; configured window {window:,} bp)."
        )

    print("[interpret] case study: computing external HARE5 feature row from hg38 interval.")
    return build_external_hare5_row(chrom38, s38, e38)


def write_hare5_case_study(
    df: pd.DataFrame,
    rank_df: pd.DataFrame,
    plot_df: pd.DataFrame,
) -> None:
    """Write HARE5 feature table and case-study distribution panel."""
    print("[interpret] building HARE5 case study panel ...")
    hare5 = find_hare5(df)
    if hare5 is None:
        print("[interpret] case study skipped (HARE5 not found or liftOver failed).")
        return
    pd.DataFrame([hare5]).to_csv("outputs/tables/hare5_features.tsv", sep="\t", index=False)

    hare5_x = transform(pd.DataFrame([hare5[FEATURE_COLS]])).iloc[0]
    top4 = rank_df.head(4)["feature"].tolist()

    fig, axes = plt.subplots(1, len(top4), figsize=(3.4 * len(top4), 3.6), sharey=False)
    for ax, feat in zip(axes, top4):
        if feat.startswith("overlaps_"):
            # bar + marker for HARE5's value
            counts = plot_df.groupby("label")[feat].mean().reindex(["CNE","HAR"])
            ax.bar(counts.index, counts.values, color=["#888","#cc4d4d"])
            ax.set_ylabel("Fraction overlapping")
            ax.axhline(hare5_x[feat], color="black", linestyle=":", lw=1)
            ax.scatter(["HAR"], [hare5_x[feat]], color="gold",
                       edgecolor="black", s=140, zorder=5,
                       label=f"HARE5 = {int(hare5_x[feat])}")
        else:
            for lab, color in [("CNE","#888"), ("HAR","#cc4d4d")]:
                vals = plot_df.loc[plot_df["label"]==lab, feat].dropna()
                ax.hist(vals, bins=40, alpha=0.55, color=color, label=lab, density=True)
            ax.axvline(hare5_x[feat], color="gold", lw=2.5, label=f"HARE5 = {hare5_x[feat]:.2f}")
        ax.set_title(PRETTY[feat], fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="best")
        sns.despine(ax=ax)
    fig.suptitle(f"{CFG['case_study'].get('alias', CFG['case_study']['name'])} within HAR/CNE distributions of top-ranked features",
                 fontsize=11, y=1.04)
    fig.tight_layout()
    fig.savefig("outputs/figures/hare5_case_study.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[interpret] -> outputs/figures/hare5_case_study.png")
    print("[interpret] -> outputs/tables/hare5_features.tsv")


def main() -> None:
    df = pd.read_csv("data/processed/features.tsv", sep="\t")

    if "--case-study-only" in sys.argv:
        rank_df = pd.read_csv("outputs/tables/top_shap_features.tsv", sep="\t")
        plot_df = transform(df[FEATURE_COLS])
        plot_df["label"] = ["HAR" if v == 1 else "CNE" for v in df["label"].values]
        write_hare5_case_study(df, rank_df, plot_df)
        return

    rf = joblib.load("outputs/models/rf.pkl")

    X_raw = df[FEATURE_COLS]
    X     = transform(X_raw)
    y     = df["label"].values

    # ---- SHAP -----------------------------------------------------------
    print("[interpret] computing SHAP values (TreeExplainer on RF) ...")
    explainer = shap.TreeExplainer(rf)
    sv = explainer.shap_values(X.values)
    # newer shap returns array of shape (n,d,2); older returns list of two arrays.
    if isinstance(sv, list):
        sv_pos = sv[1]
    elif sv.ndim == 3:
        sv_pos = sv[:, :, 1]
    else:
        sv_pos = sv  # binary classifier with single output

    # --- summary beeswarm
    plt.figure(figsize=(7.5, 5))
    shap.summary_plot(
        sv_pos, X.values, feature_names=[PRETTY[c] for c in FEATURE_COLS],
        show=False, plot_size=None,
    )
    plt.title("SHAP summary: features distinguishing HARs from matched CNEs", pad=12)
    plt.tight_layout()
    plt.savefig("outputs/figures/shap_summary.png", dpi=220, bbox_inches="tight")
    plt.close()

    # --- mean |SHAP| bar (cleaner for a 10-min talk)
    mean_abs = np.abs(sv_pos).mean(axis=0)
    order = np.argsort(mean_abs)[::-1]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(
        [PRETTY[FEATURE_COLS[i]] for i in order][::-1],
        mean_abs[order][::-1],
        color="#1f4e79",
    )
    ax.set_xlabel("Mean |SHAP value|  (impact on classifying as HAR)")
    ax.set_title("Feature importance ranking", pad=10)
    sns.despine()
    fig.tight_layout()
    fig.savefig("outputs/figures/shap_bar.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # --- ranked table with direction-of-effect via Spearman(feature, SHAP)
    rho = np.array([spearmanr(X.values[:, j], sv_pos[:, j]).statistic
                    for j in range(X.shape[1])])
    direction_label = []
    for j in order:
        r = rho[j]
        if abs(r) < 0.05:
            direction_label.append("no clear direction")
        elif r > 0:
            direction_label.append("higher value -> more HAR-like")
        else:
            direction_label.append("higher value -> more CNE-like")
    rank_df = pd.DataFrame({
        "feature":          [FEATURE_COLS[i] for i in order],
        "feature_label":    [PRETTY[FEATURE_COLS[i]] for i in order],
        "mean_abs_shap":    mean_abs[order],
        "spearman_feat_shap": rho[order],
        "direction":        direction_label,
    })

    rank_df.to_csv("outputs/tables/top_shap_features.tsv", sep="\t", index=False)
    print("[interpret] -> outputs/figures/shap_summary.png + shap_bar.png")
    print("[interpret] -> outputs/tables/top_shap_features.tsv")
    print("\nFeature ranking:\n", rank_df.to_string(index=False), "\n")

    # ---- per-feature distributions (HAR vs CNE) -------------------------
    print("[interpret] plotting feature distributions (HAR vs CNE) ...")
    plot_df = X.copy()
    plot_df["label"] = ["HAR" if v == 1 else "CNE" for v in y]
    fig, axes = plt.subplots(2, 4, figsize=(13, 6))
    for ax, col in zip(axes.flat, FEATURE_COLS):
        if col.startswith("overlaps_"):
            counts = (plot_df.groupby("label")[col]
                      .mean().reindex(["CNE","HAR"]))
            ax.bar(counts.index, counts.values,
                   color=["#888","#cc4d4d"])
            ax.set_ylabel("Fraction overlapping")
            ax.set_ylim(0, 1)
        else:
            for lab, color in [("CNE","#888"), ("HAR","#cc4d4d")]:
                vals = plot_df.loc[plot_df["label"]==lab, col].dropna()
                ax.hist(vals, bins=40, alpha=0.55, color=color, label=lab, density=True)
            ax.legend(frameon=False, fontsize=9)
        ax.set_title(PRETTY[col], fontsize=10)
        sns.despine(ax=ax)
    # hide the unused 8th panel
    axes.flat[-1].axis("off")
    fig.suptitle("Feature distributions, HARs vs matched CNEs", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig("outputs/figures/feature_distributions.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- HARE5 case study -----------------------------------------------
    write_hare5_case_study(df, rank_df, plot_df)


if __name__ == "__main__":
    sys.exit(main())
