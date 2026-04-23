"""
03_classify.py

Train two classifiers on the features table and evaluate with stratified 5-fold CV:
    1. Logistic regression (linear baseline)
    2. Random forest (interpretable nonlinear)

Outputs:
    outputs/tables/cv_metrics.tsv         per-model AUROC, AUPRC, accuracy, F1, by fold + mean
    outputs/figures/roc_curves.png        ROC + PR curves overlaid
    outputs/models/rf.pkl                 RF refit on all data (used by 04_interpret.py)
    outputs/models/lr.pkl                 LR refit on all data
    outputs/models/feature_scaler.pkl     StandardScaler used by LR

Notes:
    - Distance features are log1p-transformed before scaling (long right tails).
    - Class imbalance is handled with class_weight='balanced' (CNEs outnumber HARs ~3:1).
    - We DON'T tune hyperparameters extensively — the goal is interpretability,
      and a heavily tuned RF could overfit to spurious features.
"""
from __future__ import annotations

import sys
import yaml
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score, f1_score,
    roc_curve, precision_recall_curve,
)
from pathlib import Path

CFG = yaml.safe_load(open("config.yaml"))
SEED = CFG["modeling"]["random_seed"]
N_SPLITS = CFG["modeling"]["cv_folds"]

FEATURE_COLS = [
    "gc_content", "length", "mean_phastcons",
    "dist_to_nearest_tss", "dist_to_nearest_brain_tss",
    "overlaps_ccre", "overlaps_fetal_brain_enh",
]
LOG_FEATURES = ["length", "dist_to_nearest_tss", "dist_to_nearest_brain_tss"]


def transform_features(X: pd.DataFrame) -> pd.DataFrame:
    """Log1p-transform skewed distance/length features. Returns a new frame."""
    Xt = X.copy()
    for c in LOG_FEATURES:
        Xt[c] = np.log1p(Xt[c])
    return Xt


def main() -> None:
    df = pd.read_csv("data/processed/features.tsv", sep="\t")
    X_raw = df[FEATURE_COLS]
    X = transform_features(X_raw)
    y = df["label"].values

    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    rows = []

    # for plotting mean ROC/PR
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for model_name, build_model in [
        ("logreg", lambda: LogisticRegression(
            max_iter=5000, class_weight="balanced", random_state=SEED)),
        ("rf",     lambda: RandomForestClassifier(
            n_estimators=CFG["modeling"]["rf_n_estimators"],
            max_depth=CFG["modeling"]["rf_max_depth"],
            class_weight="balanced",
            random_state=SEED, n_jobs=-1)),
    ]:
        fold_tprs = []
        mean_fpr = np.linspace(0, 1, 200)
        fold_precs = []
        mean_rec = np.linspace(0, 1, 200)

        for fold, (tr, te) in enumerate(skf.split(X, y), start=1):
            scaler = StandardScaler().fit(X.iloc[tr])
            Xtr = scaler.transform(X.iloc[tr]) if model_name == "logreg" else X.iloc[tr].values
            Xte = scaler.transform(X.iloc[te]) if model_name == "logreg" else X.iloc[te].values

            m = build_model().fit(Xtr, y[tr])
            p = m.predict_proba(Xte)[:, 1]
            yhat = (p >= 0.5).astype(int)
            rows.append({
                "model": model_name, "fold": fold,
                "auroc":    roc_auc_score(y[te], p),
                "auprc":    average_precision_score(y[te], p),
                "accuracy": accuracy_score(y[te], yhat),
                "f1":       f1_score(y[te], yhat),
            })
            # ROC
            fpr, tpr, _ = roc_curve(y[te], p)
            fold_tprs.append(np.interp(mean_fpr, fpr, tpr))
            fold_tprs[-1][0] = 0.0
            # PR
            prec, rec, _ = precision_recall_curve(y[te], p)
            fold_precs.append(np.interp(mean_rec, rec[::-1], prec[::-1]))

        # plot mean curves
        mean_tpr = np.mean(fold_tprs, axis=0); mean_tpr[-1] = 1.0
        mean_prec = np.mean(fold_precs, axis=0)
        axes[0].plot(mean_fpr, mean_tpr,
                     label=f"{model_name} (AUROC={np.mean([r['auroc'] for r in rows if r['model']==model_name]):.3f})",
                     lw=2)
        axes[1].plot(mean_rec, mean_prec,
                     label=f"{model_name} (AUPRC={np.mean([r['auprc'] for r in rows if r['model']==model_name]):.3f})",
                     lw=2)

    # finalize plot
    axes[0].plot([0, 1], [0, 1], "k--", alpha=0.4)
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC (5-fold CV mean)")
    axes[0].legend(loc="lower right", frameon=False)

    pos_rate = y.mean()
    axes[1].axhline(pos_rate, color="k", linestyle="--", alpha=0.4,
                    label=f"baseline ({pos_rate:.2f})")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision–Recall (5-fold CV mean)")
    axes[1].legend(loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig("outputs/figures/roc_curves.png", dpi=200)
    plt.close(fig)
    print("[classify] -> outputs/figures/roc_curves.png")

    # save metrics
    metrics = pd.DataFrame(rows)
    metrics_summary = (metrics.groupby("model")
                       [["auroc","auprc","accuracy","f1"]]
                       .agg(["mean","std"]).round(3))
    metrics.to_csv("outputs/tables/cv_metrics.tsv", sep="\t", index=False)
    metrics_summary.to_csv("outputs/tables/cv_metrics_summary.tsv", sep="\t")
    print("[classify] -> outputs/tables/cv_metrics{,_summary}.tsv")
    print("\n", metrics_summary, "\n")

    # refit on all data and persist (used by 04_interpret.py)
    Path("outputs/models").mkdir(parents=True, exist_ok=True)
    scaler = StandardScaler().fit(X)
    joblib.dump(scaler, "outputs/models/feature_scaler.pkl")

    lr = LogisticRegression(max_iter=5000, class_weight="balanced", random_state=SEED)
    lr.fit(scaler.transform(X), y)
    joblib.dump(lr, "outputs/models/lr.pkl")

    rf = RandomForestClassifier(
        n_estimators=CFG["modeling"]["rf_n_estimators"],
        max_depth=CFG["modeling"]["rf_max_depth"],
        class_weight="balanced",
        random_state=SEED, n_jobs=-1,
    )
    rf.fit(X.values, y)
    joblib.dump(rf, "outputs/models/rf.pkl")
    print("[classify] models -> outputs/models/{lr,rf,feature_scaler}.pkl")


if __name__ == "__main__":
    sys.exit(main())
