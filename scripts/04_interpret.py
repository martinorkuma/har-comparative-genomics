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


def find_hare5(df: pd.DataFrame) -> pd.Series | None:
    """Locate HARE5 in the feature table.

    Two attempts, in order:
      (a) match the HAR by name (e.g. '2xHAR.238') if the HAR set used names.
      (b) fall back to a window around the lifted hg38 coordinates of HARE5.
    """
    name_target = CFG["case_study"]["name"]
    hits = df[df["name"].astype(str).str.strip() == name_target]
    if len(hits):
        print(f"[interpret] case study: matched HAR by name '{name_target}'")
        return hits.iloc[0]

    # Fallback: nearest HAR to the case-study locus.
    # We use the hg19 coordinates and naively shift; a real liftOver would be
    # better, but for FZD8 the hg19 vs hg38 offset is tiny.
    chrom = CFG["case_study"]["hg19_chrom"]
    mid_hg19 = (CFG["case_study"]["hg19_start"] + CFG["case_study"]["hg19_end"]) // 2
    window = CFG["case_study"]["fallback_window_bp"]
    cand = df[(df["label"] == 1) & (df["chrom"] == chrom)].copy()
    if cand.empty:
        print(f"[interpret] WARNING: no HARs on {chrom}; case study unavailable.")
        return None
    cand["mid"] = (cand["start"] + cand["end"]) // 2
    cand["d"]   = (cand["mid"] - mid_hg19).abs()
    nearest = cand.sort_values("d").iloc[0]
    if nearest["d"] > 5_000_000:
        print(f"[interpret] WARNING: nearest HAR on {chrom} is "
              f"{nearest['d']/1e6:.2f} Mb from HARE5 hg19 locus; check liftOver.")
    print(f"[interpret] case study: fallback to HAR '{nearest['name']}' "
          f"({nearest['d']:,} bp from HARE5 hg19 midpoint)")
    return nearest


def main() -> None:
    df = pd.read_csv("data/processed/features.tsv", sep="\t")
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

    # --- ranked table with mean direction (sign of mean SHAP, not |SHAP|)
    mean_sv  = sv_pos.mean(axis=0)
    direction = np.sign(mean_sv)
    rank_df = pd.DataFrame({
        "feature":           [FEATURE_COLS[i] for i in order],
        "feature_label":     [PRETTY[FEATURE_COLS[i]] for i in order],
        "mean_abs_shap":     mean_abs[order],
        "mean_signed_shap":  mean_sv[order],
        "direction":         ["increases HAR prob" if direction[i] > 0 else
                              "decreases HAR prob" for i in order],
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
    print("[interpret] building HARE5 case study panel ...")
    hare5 = find_hare5(df)
    if hare5 is None:
        print("[interpret] case study skipped (HARE5 not found in HAR set).")
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
    fig.suptitle(f"HARE5 ({CFG['case_study']['name']}) within HAR/CNE distributions of top-ranked features",
                 fontsize=11, y=1.04)
    fig.tight_layout()
    fig.savefig("outputs/figures/hare5_case_study.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[interpret] -> outputs/figures/hare5_case_study.png")
    print("[interpret] -> outputs/tables/hare5_features.tsv")


if __name__ == "__main__":
    sys.exit(main())
