# Poster Outline

**Working title:** *What makes a HAR a HAR? Interpretable machine learning on human-accelerated noncoding regions*

**Poster format:** 48" × 36" landscape, three-column layout.

**Central message:** Human accelerated regions (HARs) can be compared against carefully matched conserved non-accelerated elements (CNEs) using interpretable machine learning. The poster should make the SHAP result the visual anchor, then use distributions, sensitivity checks, and the HARE5/FZD8 case study to explain and qualify that result.

---

## Visual hierarchy

The reader's eye should land first on the SHAP beeswarm plot. The left column sets up the biological question and computational workflow. The middle column presents the main model interpretation and supporting distributions. The right column translates the model result into a biologically recognizable case study and concise takeaways.

```text
┌─────────────────────┬─────────────────────────────┬─────────────────────┐
│  1. Background      │   3. SHAP feature ranking   │  5. HARE5 / FZD8    │
│     & question      │      CENTERPIECE            │     case study      │
│                     │                             │                     │
├─────────────────────┤                             ├─────────────────────┤
│  2. Approach        │   4. Feature distributions  │  6. Conclusions &   │
│     workflow +      │      HAR vs CNE             │     future work     │
│     features        │                             │                     │
│                     ├─────────────────────────────┤                     │
│                     │  4b. Sensitivity: brain-TSS │                     │
│                     │      proximity, 3-way       │                     │
│                     │  4c. Classifier performance │                     │
│                     │      AUROC/AUPRC/F1 strip   │                     │
└─────────────────────┴─────────────────────────────┴─────────────────────┘
```

### Recommended poster flow

1. **Panel 1:** What are HARs, and why ask this question?
2. **Panel 2:** How were HARs compared to matched CNEs?
3. **Panel 3:** Which features distinguish HARs from controls?
4. **Panel 4 / 4b / 4c:** Do raw distributions, sensitivity checks, and classifier metrics support the interpretation?
5. **Panel 5:** How does the textbook HARE5/FZD8 example fit into the feature space?
6. **Panel 6:** What is the biological interpretation, and what should be done next?

---

## Figure inventory

| Poster panel | Graphic | Repo target path | Purpose |
|---|---|---|---|
| 1 | HAR concept schematic | `outputs/figures/poster_har_concept_schematic.png` | Explain conservation across non-human species followed by human-lineage acceleration |
| 2 | Workflow diagram | `outputs/figures/poster_workflow_diagram.png` | Show data → features → models → SHAP interpretation |
| 2 | Feature rationale mini-table | Built into poster, or exported as `outputs/figures/poster_feature_rationale_table.png` | Connect each feature to a biological rationale |
| 3 | SHAP beeswarm | `outputs/figures/shap_summary.png` | Main result; visual centerpiece |
| 3 | SHAP bar inset | `outputs/figures/shap_bar.png` | Quick feature ranking |
| 4 | Feature distributions | `outputs/figures/feature_distributions.png` | Show whether raw HAR/CNE distributions agree with SHAP |
| 4b | Brain-TSS sensitivity violin | `outputs/figures/sensitivity_brain_tss.png` | Test whether brain-gene proximity result is robust or control-driven |
| 4c | ROC/PR curves | `outputs/figures/roc_curves.png` | Show predictive signal without overemphasizing model accuracy |
| 4c | CV metrics table | `outputs/tables/cv_metrics_summary.tsv` → poster table | Summarize AUROC, AUPRC, and F1 |
| 5 | HARE5/FZD8 case study | `outputs/figures/hare5_case_study.png` | Place external HARE5 interval in HAR/CNE feature context |
| 6 | Summary model graphic | `outputs/figures/poster_takeaway_model.png` | Simple visual synthesis of final interpretation |

---

## Panel 1 — Background & question

### Poster text

**Human accelerated regions (HARs)** are short genomic elements that are deeply conserved across non-human vertebrates but show unusually rapid sequence change on the human lineage. Many HARs are noncoding and are hypothesized to act as developmental regulatory elements, especially in contexts relevant to neurodevelopment.

**Question:** Which genomic and regulatory features best distinguish HARs from matched conserved non-accelerated elements?

**Hypothesis:** HARs are better distinguished by regulatory context than by simple sequence features alone.

### Graphic 1A — HAR concept schematic

**Target file:** `outputs/figures/poster_har_concept_schematic.png`

**Design:** A small alignment-style schematic with four rows:

```text
Mouse      A C T G C A T G C T A C
Dog        A C T G C A T G C T A C
Chimp      A C T G C A T G C T A C
Human      A T T A C G T T C A A C
                 └── HAR ──┘
```

Add a label above the highlighted block:

```text
Conserved across non-human species → accelerated on human lineage
```

**Visual encoding:**

- Conserved non-human sequence: neutral gray
- Human substitutions: warm red
- HAR interval box: warm red outline
- Optional genome track below: small “noncoding regulatory element” label

**Caption:**  
*HARs are conserved across non-human species but show accelerated sequence divergence on the human lineage, motivating tests of whether this divergence occurs in distinctive regulatory contexts.*

**References to include in this panel:** Pollard 2006; Capra 2013; Doan 2016; Keough 2023.

---

## Panel 2 — Approach

### Poster text

This project compares published HARs to matched CNEs, engineers a focused set of genomic and regulatory features, trains interpretable classifiers, and uses SHAP to identify which features most strongly separate HARs from controls.

### Graphic 2A — Workflow diagram

**Target file:** `outputs/figures/poster_workflow_diagram.png`

**Recommended structure:**

```text
Published HARs
     │
     ├── Match to CNEs
     │      length ±25%
     │      mean phastCons ±0.05
     │      3 CNEs per HAR
     │      ±10 kb HAR exclusion buffer
     │
     ▼
Feature engineering
7 genomic + regulatory features
     │
     ▼
Models
Logistic regression + random forest
5-fold cross-validation
     │
     ▼
Interpretation
SHAP ranking + raw distributions + sensitivity checks
```

**Design notes:**

- Use a horizontal or vertical pipeline depending on available space.
- Keep each step as a rounded rectangle.
- Use arrows thick enough to remain readable at print scale.
- Add a small “matched-control design” badge near the CNE matching step.

**Caption:**  
*Matched CNEs provide a conservative background for testing whether HARs occupy distinctive regulatory contexts rather than merely reflecting length or conservation differences.*

### Graphic 2B — Feature rationale box

**Target:** Build directly into poster as a compact table.

| Feature | Biological rationale |
|---|---|
| Element length | Controls for interval size and annotation opportunity |
| Mean phastCons | Controls for baseline evolutionary conservation |
| GC content | Captures local sequence composition |
| Regulatory overlap | Tests enrichment in annotated enhancer/promoter-like regions |
| Distance to nearest TSS | Measures broad gene-proximity context |
| Distance to nearest brain-expressed gene/TSS | Tests the neurodevelopmental proximity hypothesis |
| Local transposable element density | Captures repeat-rich regulatory neighborhood context |

**Caption:**  
*The feature set intentionally mixes control variables, sequence composition, gene proximity, and regulatory-context annotations.*

### QC note

Add a small QC box using values from:

```text
outputs/tables/cne_matching_qc.tsv
```

Recommended text template:

```text
Matching QC: HARs and CNEs were similar in length and mean phastCons after matching, supporting interpretation of downstream regulatory-context differences.
```

Replace with exact numbers once available:

```text
n HARs = ___; n CNEs = ___
median length HAR/CNE = ___ / ___ bp
median phastCons HAR/CNE = ___ / ___
```

---

## Panel 3 — SHAP feature ranking

**Role:** Centerpiece. This should be the largest figure on the poster.

### Main graphic — SHAP beeswarm

**Target file:** `outputs/figures/shap_summary.png`

Use the beeswarm, not only the bar plot. It shows both feature importance and directionality.

**Placement:** Top-middle, occupying most of the center column.

**Caption template:**  
*SHAP interpretation identifies [TOP FEATURE 1] and [TOP FEATURE 2] as the strongest contributors to HAR/CNE classification, indicating that HARs are distinguished primarily by [BIOLOGICAL INTERPRETATION], not by classification accuracy alone.*

Replace bracketed values after reviewing final SHAP output.

### Inset — SHAP bar plot

**Target file:** `outputs/figures/shap_bar.png`

Place as a small inset in the upper-right corner of the SHAP panel.

**Caption:**  
*Bar summary confirms the same top-ranked features in a simpler at-a-glance format.*

### Annotation arrows

Add one or two arrows directly on the beeswarm:

```text
Top feature: strongest model signal
Direction: higher values push prediction toward HAR or CNE
```

Use the model’s actual SHAP direction:

```text
High [feature] → HAR
Low [feature] → CNE
```

or

```text
High [feature] → CNE
Low [feature] → HAR
```

Do not write the direction until the final beeswarm is inspected.

---

## Panel 4 — Feature distributions

### Main graphic

**Target file:** `outputs/figures/feature_distributions.png`

Use compact HAR vs CNE histograms, density plots, or bar plots for the features most relevant to the SHAP result.

**Purpose:** This panel should let readers verify whether the model’s interpretation is visible in the raw data.

**Caption template:**  
*Raw feature distributions support the SHAP ranking: [TOP FEATURE] differs visibly between HARs and matched CNEs, while matched variables such as length and conservation show reduced separation.*

### Design notes

- Do not show every feature if space is tight.
- Prioritize:
  1. Top SHAP feature
  2. Second SHAP feature
  3. Length
  4. Mean phastCons
  5. Brain-TSS distance or regulatory overlap, depending on model result
- Use the same colors throughout:
  - HAR: `#cc4d4d`
  - CNE: `#888888`

---

## Panel 4b — Sensitivity: brain-TSS proximity

### Graphic

**Target file:** `outputs/figures/sensitivity_brain_tss.png`

Use a three-way violin or box/violin hybrid:

```text
HAR vs matched CNE vs random length-matched control
```

**Purpose:** This panel should determine whether the brain-TSS proximity signal reflects HAR biology, the CNE matching design, or anomalously proximal CNEs.

### Required table reference

Use:

```text
outputs/tables/sensitivity_brain_tss.tsv
```

Report medians and Mann-Whitney U p-values inline.

### Caption templates

Choose the correct caption after inspecting the sensitivity output:

**Matching artifact interpretation:**  
*Brain-TSS proximity appears sensitive to the matched-control design: matched CNEs are unusually positioned relative to brain-expressed genes compared with random length-matched controls.*

**Real reversal interpretation:**  
*HARs are farther from brain-expressed TSSs than matched CNEs, suggesting that linear proximity alone does not capture HAR regulatory targeting.*

**Anomalously proximal CNE interpretation:**  
*Matched CNEs are unusually close to brain-expressed TSSs, indicating that the apparent HAR/CNE difference is driven partly by the control set.*

---

## Panel 4c — Classifier performance strip

### Graphic

**Target files:**

```text
outputs/figures/roc_curves.png
outputs/tables/cv_metrics_summary.tsv
```

Place ROC and PR curves as a thin strip below Panel 4b. Add a tiny metrics table beside or beneath the curves.

### Metrics table format

| Model | AUROC | AUPRC | F1 |
|---|---:|---:|---:|
| Logistic regression | mean ± SD | mean ± SD | mean ± SD |
| Random forest | mean ± SD | mean ± SD | mean ± SD |

### Caption

*Both models achieve modest but above-chance discrimination with consistent feature rankings; the project emphasizes biological interpretation over maximizing predictive accuracy.*

### Design note

Keep this section visually smaller than SHAP. The performance result supports credibility but should not become the main story.

---

## Panel 5 — Case study: HARE5 / FZD8

### Poster text

HARE5 is a well-known human accelerated enhancer near **FZD8**, a gene implicated in brain development. Because HARE5/2xHAR.238 is not included in the modeled Doan HAR table, the lifted-over Boyd interval is scored as an external reference using the same feature sources rather than treated as a training-set observation.

### Graphic

**Target file:** `outputs/figures/hare5_case_study.png`

**Recommended design:** A compact multi-part graphic:

```text
A. Genome context
   HARE5 interval ───── nearby FZD8

B. Feature-space placement
   HARE5 marker over HAR/CNE distributions

C. External-reference note
   Scored with same features, excluded from model training
```

**Visual encoding:**

- HAR distribution: warm red
- CNE distribution: gray
- HARE5 marker: gold `#d4a017`
- FZD8 gene label: black or dark gray

**Caption:**  
*HARE5/FZD8 provides an external biological reference: the interval is scored using the same feature definitions and placed into the HAR/CNE feature space without being used to train the classifier.*

### Design note

Do not overclaim from a single case. Present HARE5 as an interpretive anchor, not as proof of the full model.

---

## Panel 6 — Conclusions & future work

### Main conclusion bullets

1. **Result:** The strongest SHAP-ranked feature(s) distinguishing HARs from matched CNEs are **[TOP FEATURE 1]** and **[TOP FEATURE 2]**.
2. **Biological interpretation:** HARs preferentially show **[REGULATORY CONTEXT INTERPRETATION]**, while linear distance to brain-expressed genes does not support a simple “near brain genes” model under this matched-control design.
3. **Sensitivity check:** The three-way brain-TSS analysis clarifies whether the observed proximity pattern reflects HAR biology, the conservation-matched CNE background, or an artifact of control selection.
4. **Limitations:** Remaining issues include residual conservation signal, regulatory annotation bias, and the fact that nearest-gene distance is an imperfect proxy for enhancer target genes.
5. **Next steps:** Extend the feature set with 3D chromatin contacts, developmental brain regulatory atlases, and single-cell brain enhancer annotations.

### Graphic 6A — Summary model graphic

**Target file:** `outputs/figures/poster_takeaway_model.png`

**Recommended structure:**

```text
HAR classification signal
        │
        ├── Regulatory annotation context
        ├── Conservation-matched comparison
        ├── Brain-gene proximity sensitivity
        └── HARE5 external reference
```

**Caption:**  
*The strongest interpretation comes from combining model explanation, raw feature distributions, sensitivity analysis, and an external biological reference case.*

---

## Color and style guide

### Core colors

| Element | Color |
|---|---|
| HAR | Warm red `#cc4d4d` |
| CNE | Neutral gray `#888888` |
| HARE5 | Gold `#d4a017` |
| Text | Near-black `#222222` |
| Panel background | White or very light gray |
| Gridlines | Light gray, low contrast |

### Typography

- Use one sans-serif typeface throughout: Inter, Source Sans, Arial, or department template font.
- Avoid more than two heading levels.
- Suggested minimum sizes for 48" × 36":
  - Title: 90–110 pt
  - Panel headers: 44–56 pt
  - Body text: 28–34 pt
  - Figure captions: 24–28 pt
  - Axis labels: ≥24 pt
  - Tick labels: ≥20–24 pt, preferably ≥24 pt

### Figure rules

- Every figure should have a caption that states the biological punchline.
- Avoid captions that only describe the axes.
- Use direct labels where possible instead of legends.
- Keep HAR/CNE/HARE5 colors consistent across all graphics.
- Do not use gold except for HARE5.

---

## Suggested scripts to add

These helper scripts would make the poster graphics reproducible.

```text
scripts/make_poster_figures.sh
src/visualization/make_poster_schematics.py
src/visualization/make_poster_tables.py
```

### `scripts/make_poster_figures.sh`

Purpose: one-command poster graphic regeneration.

```bash
#!/usr/bin/env bash
set -euo pipefail

python src/visualization/make_poster_schematics.py
python src/visualization/make_poster_tables.py

echo "Poster graphics written to outputs/figures/"
```

### `src/visualization/make_poster_schematics.py`

Purpose: create schematic-only figures that are not direct model outputs.

Recommended outputs:

```text
outputs/figures/poster_har_concept_schematic.png
outputs/figures/poster_workflow_diagram.png
outputs/figures/poster_takeaway_model.png
```

### `src/visualization/make_poster_tables.py`

Purpose: convert TSV outputs into small poster-ready tables.

Recommended inputs:

```text
outputs/tables/cne_matching_qc.tsv
outputs/tables/cv_metrics_summary.tsv
outputs/tables/sensitivity_brain_tss.tsv
```

Recommended outputs:

```text
outputs/figures/poster_matching_qc_table.png
outputs/figures/poster_cv_metrics_table.png
outputs/figures/poster_brain_tss_sensitivity_table.png
```

---

## References to include

Keep references compact. Use 5–7 key references maximum on the poster.

Recommended core references:

- Pollard et al. 2006 — original HAR framing
- Capra et al. 2013 — HARs and developmental enhancers
- Doan et al. 2016 — HAR mutations, cognition, social behavior
- Boyd et al. 2015 — HARE5/FZD8 enhancer case study
- Keough et al. 2023 — 3D genome context of HAR loci
- Cui et al. 2025 — recent neuron-focused HAR characterization, if space allows

---

## Pre-print checklist

- [ ] Every figure is at least 300 dpi.
- [ ] SHAP beeswarm is the largest figure on the poster.
- [ ] Every figure has a caption stating the punchline.
- [ ] Every claim that requires a citation has one.
- [ ] HARE5 appears by name in Panel 1, Panel 5, and Panel 6.
- [ ] Panel 4b caption explicitly names the supported interpretation.
- [ ] No body text is smaller than 24 pt.
- [ ] HAR, CNE, and HARE5 colors are consistent across all figures.
- [ ] ROC/PR performance is visually subordinate to SHAP interpretation.
- [ ] Print at 25% scale on tabloid paper to sanity-check legibility.
- [ ] Confirm that no model-training figure treats HARE5 as a training-set observation.
- [ ] Confirm that all output paths in the poster match the actual repo paths.

---

## Final poster assembly checklist

Before exporting the final PDF:

```bash
# From repo root
ls outputs/figures/shap_summary.png
ls outputs/figures/shap_bar.png
ls outputs/figures/feature_distributions.png
ls outputs/figures/sensitivity_brain_tss.png
ls outputs/figures/roc_curves.png
ls outputs/figures/hare5_case_study.png
ls outputs/tables/cne_matching_qc.tsv
ls outputs/tables/cv_metrics_summary.tsv
ls outputs/tables/sensitivity_brain_tss.tsv
```

Then verify:

```bash
git status --short
```

Expected tracked update:

```text
M docs/poster_outline.md
```

Expected untracked/generated outputs, if created locally:

```text
outputs/figures/...
outputs/tables/...
```
