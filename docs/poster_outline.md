# Poster Outline

**Title:** What makes a HAR a HAR? Interpretable ML on human-accelerated noncoding regions

**Format:** Standard 48"×36" landscape (or institution-required dimensions). Three columns.

---

## Visual hierarchy

The reader's eye should land on the SHAP plot. Everything else either sets it up
(left column) or supports / qualifies it (right column).

```
┌─────────────────────┬─────────────────────────────┬─────────────────────┐
│  1. Background      │   3. SHAP feature ranking   │  5. HARE5 / FZD8    │
│     & question      │      (CENTERPIECE)          │     case study      │
│                     │                             │                     │
├─────────────────────┤                             ├─────────────────────┤
│  2. Approach        │   4. Feature distributions  │  6. Conclusions &   │
│     (workflow +     │      (HAR vs CNE, small)    │     future work     │
│      features)      │                             │                     │
│                     ├─────────────────────────────┤                     │
│                     │  4b. Sensitivity: brain-TSS │                     │
│                     │      proximity (3-way)      │                     │
│                     │  4c. Classifier perf strip  │                     │
│                     │      (AUROC/AUPRC table)    │                     │
└─────────────────────┴─────────────────────────────┴─────────────────────┘
```

---

## Panel 1 — Background & question (top-left)

- Two short paragraphs:
  - What HARs are, why they matter for human evolution.
  - The open question this project addresses.
- One small schematic: vertebrate alignment with a HAR highlighted, showing
  conservation across non-humans then human-lineage divergence.
- Key references inline: Pollard 2006, Capra 2013, Doan 2016.

## Panel 2 — Approach (bottom-left)

- Pipeline diagram: HAR set + matched CNEs → 7 features → LR + RF (5-fold CV) → SHAP.
- Box listing the 7 features with one-line biological rationale each.
- Note on matching strategy: length (±25%) and conservation (±0.05 mean phastCons),
  3:1 CNE:HAR ratio, ±10 kb buffer to keep CNEs distant from HARs.
- QC bullet: numbers from `outputs/tables/cne_matching_qc.tsv` (n, length and
  conservation summaries showing the match worked).

## Panel 3 — SHAP feature ranking (CENTERPIECE, top-middle, large)

- Use the **beeswarm** here (`outputs/figures/shap_summary.png`), not the bar.
  Beeswarm is information-dense and rewards the close inspection a poster invites
  but a talk audience doesn't have time for.
- Inset (small): the bar version (`outputs/figures/shap_bar.png`) for at-a-glance
  ranking.
- One-sentence caption that tells the reader the punchline directly.
- Annotation arrows on the top 1-2 features pointing to a short interpretive note.

## Panel 4 — Feature distributions (middle column, below SHAP)

- Use `outputs/figures/feature_distributions.png` (HAR vs CNE histograms / bars).
- Highlights *why* the SHAP ranking comes out the way it does — readers can verify
  the model's claims against the raw distributions.
- Small caption only; the SHAP panel does the talking.

## Panel 4b — Sensitivity: brain-TSS proximity (small box, middle column)

- Three-way violin from `outputs/figures/sensitivity_brain_tss.png`: HAR (red),
  matched CNE (gray), random length-matched (blue).
- One-sentence caption stating which pre-decided interpretation the data
  support (matching artifact / real reversal / anomalously proximal CNEs).
- Cite `outputs/tables/sensitivity_brain_tss.tsv` for medians and Mann-Whitney
  U p-values inline.

## Panel 4c — Classifier performance strip (smaller, below 4b)

- ROC + PR curves (`outputs/figures/roc_curves.png`), single row.
- Tiny table: AUROC, AUPRC, F1 for LR vs RF (mean ± SD across 5 folds) from
  `outputs/tables/cv_metrics_summary.tsv`.
- One sentence: *"Both models achieve modest but well-above-chance discrimination
  with consistent feature rankings; the project emphasizes interpretation over
  predictive accuracy."*

## Panel 5 — Case study: HARE5 / FZD8 (top-right)

- Motivation: Boyd 2015's HARE5/FZD8 work is the textbook example of why HARs matter.
- HARE5/2xHAR.238 is not in the modeled Doan HAR table, so `04_interpret.py` scores the 
  lifted-over Boyd interval as an external reference using the same feature sources.
- Show `hare5_case_study.png`. Punchline: place the real HARE5 interval within the 
  HAR/CNE feature distributions without treating it as a training-set observation.

## Panel 6 — Conclusions & future (bottom-right)

- Three bullets, in this order:
  1. **Result:** the top SHAP feature(s) most distinguishing HARs from matched CNEs.
  2. **Biological interpretation:** HARs preferentially overlap annotated regulatory
     elements, while linear TSS distance does not support a simple "near brain genes"
     story in this matched-control design — and the three-way sensitivity check
     (Panel 4b) clarifies whether this reflects HAR biology or the conservation
     matching itself. 
  3. **Limitations & next steps:** residual phastCons signal, annotation bias,
     distance-to-gene as a proxy for
     regulatory target, opportunity to extend with 3D contacts (Keough 2023) or
     single-cell brain regulatory atlases.
- Acknowledgments + references at the bottom.

---

## Color and style

- HAR vs CNE color scheme: **HAR = warm red (#cc4d4d)**, **CNE = neutral gray (#888)**.
  Used consistently across every figure.
- HARE5 highlight color: **gold (#d4a017)**. Used only for HARE5, never for anything
  else.
- One sans-serif typeface throughout (Inter, Source Sans, or your dept template).
- Avoid more than two heading levels.

## Pre-print checklist

- [ ] Every figure is at least 300 dpi.
- [ ] Every figure has a caption that states the punchline, not just labels axes.
- [ ] Every claim that requires a citation has one.
- [ ] HARE5 appears by name in at least Panel 1 (motivation), Panel 5 (case study),
      and Panel 6 (conclusion).
- [ ] Sensitivity panel (4b) caption names the supported interpretation explicitly,
      not just "see table".
- [ ] No panel uses font smaller than 24 pt for body text.
- [ ] Print at 25% scale on tabloid paper to sanity-check legibility.
