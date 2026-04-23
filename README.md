# HAR Comparative Genomics

**Do Human Accelerated Regions (HARs) sit in distinctive genomic neighborhoods compared to
matched conserved non-accelerated elements (CNEs)? An interpretable-ML approach.**

BCG 540 — Functional Genomics — Spring 2026

---

## What this project does (in one sentence)

Train a small interpretable classifier on hand-engineered genomic features for HARs vs.
length- and conservation-matched CNEs, then use SHAP to identify *which* features
distinguish HARs — and verify the top-ranked feature with a worked example at the
**HARE5 / FZD8** locus (Boyd et al. 2015).

## Scope discipline

This project was deliberately tightened after instructor feedback that an end-to-end
data → model → interpretation pipeline is too much for a 10-minute talk. The pipeline
still runs end-to-end (the poster and paper need it), but the *deliverables* are organized
around **one central result**:

> A SHAP-ranked list of features that distinguish HARs from matched CNEs, anchored by
> the HARE5/FZD8 case study.

Everything else — classifier metrics, ROC curves, feature distributions — is supporting
material on the poster, not the talk.

## Pipeline (one command)

```bash
bash setup.sh                  # one-time: conda env + tool checks
conda activate har-ml
bash scripts/run_all.sh        # end-to-end pipeline
```

`run_all.sh` runs four stages in order:

| Stage | Script | What it does | Key output |
|---|---|---|---|
| 1 | `01_acquire.sh` | Download HARs, phastCons CNEs, ENCODE cCREs, GENCODE, GTEx brain expression. LiftOver HARs to hg38 if needed. | `data/raw/`, `data/processed/hars.hg38.bed`, `data/processed/cnes.hg38.bed` |
| 2 | `02_build_features.py` | Compute 7 features per element (HAR or CNE) → one tidy table. | `data/processed/features.tsv` |
| 3 | `03_classify.py` | Logistic regression (baseline) + Random Forest with stratified 5-fold CV. | `outputs/tables/cv_metrics.tsv`, `outputs/models/rf.pkl`, ROC figure |
| 4 | `04_interpret.py` | SHAP values on RF; case-study panel for HARE5. | `outputs/figures/shap_summary.png`, `outputs/figures/hare5_case_study.png` |

## Features (deliberately small set)

Seven features per element. Each is justified by a specific biological hypothesis about
*why* it might differ between HARs and CNEs.

| # | Feature | Hypothesis if HAR-distinguishing |
|---|---|---|
| 1 | GC content | HARs may show compositional bias from accelerated substitution. |
| 2 | Element length | Length affects mutation opportunity; matched in controls but kept as a sanity check. |
| 3 | phastCons 100-way score | HARs are by definition deeply conserved before acceleration; should match CNEs (control). |
| 4 | Distance to nearest TSS (any GENCODE gene) | Tests whether HARs sit near genes generally. |
| 5 | **Distance to nearest brain-expressed TSS** (GTEx brain median TPM > 5) | **Core hypothesis test:** HARs preferentially flank brain-relevant genes. |
| 6 | Overlap with ENCODE cCRE (any) | Tests whether HARs are more likely to be annotated regulatory. |
| 7 | Overlap with fetal-brain-active enhancer (Roadmap E081/E082) | More targeted regulatory test. |

## The concrete example: HARE5 → FZD8

`04_interpret.py` produces a dedicated case-study figure for **HARE5 (2xHAR.238)**, the
HAR shown in Boyd et al. (2015, *Current Biology*) to drive earlier and broader cortical
progenitor expression of *FZD8* and produce a measurably larger neocortex in transgenic
mice. The figure shows where HARE5's feature values sit within the full HAR and CNE
distributions — making concrete what the SHAP plot says in the abstract.

If your top SHAP feature ends up *not* being one HARE5 illustrates well, swap the case
study by editing `CASE_STUDY_HAR` in `config.yaml`. Other good options: HAR1F (Pollard
2006, Cajal–Retzius neurons), 2xHAR.114 (cortical enhancer activity).

## Repository layout

```
har-comparative-genomics/
├── README.md
├── setup.sh                    # one-time env setup
├── environment.yml             # conda
├── requirements.txt            # pip extras
├── config.yaml                 # all paths, URLs, parameters
├── .gitignore
├── scripts/
│   ├── 01_acquire.sh           # data download + liftOver + matched CNE construction
│   ├── 02_build_features.py    # feature engineering
│   ├── 03_classify.py          # LR + RF + 5-fold CV
│   ├── 04_interpret.py         # SHAP + HARE5 case study
│   └── run_all.sh              # orchestrator
├── data/
│   ├── raw/                    # downloaded files (gitignored)
│   └── processed/              # derived (gitignored)
├── outputs/
│   ├── figures/                # poster figures
│   ├── tables/                 # poster tables
│   └── models/                 # trained classifiers
├── logs/                       # per-stage logs (gitignored)
└── docs/
    ├── poster_outline.md       # poster panel-by-panel
    ├── talk_outline.md         # 10-min talk slide-by-slide
    └── paper_outline.md        # written paper structure
```

## Deliverables

- **Poster** (primary). See `docs/poster_outline.md`.
- **10-minute talk.** See `docs/talk_outline.md`.
- **Written paper.** See `docs/paper_outline.md`.

## References

- Boyd, J. L. et al. (2015). Human-chimpanzee differences in a FZD8 enhancer alter cell-cycle dynamics in the developing neocortex. *Current Biology* 25:772–779.
- Capra, J. A., Erwin, G. D., McKinsey, G., Rubenstein, J. L. R., & Pollard, K. S. (2013). Many human accelerated regions are developmental enhancers. *Phil. Trans. R. Soc. B* 368:20130025.
- Doan, R. N. et al. (2016). Mutations in human accelerated regions disrupt cognition and social behavior. *Cell* 167:341–354.
- Keough, K. C. et al. (2023). Three-dimensional genome rewiring in loci with human accelerated regions. *Science* 380:eabm1696.
- Pollard, K. S. et al. (2006). An RNA gene expressed during cortical development evolved rapidly in humans. *Nature* 443:167–172.
