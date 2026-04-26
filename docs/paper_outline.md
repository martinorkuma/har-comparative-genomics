# Written Paper Outline

**Target:** ~6-8 pages, single column, 11 pt, 1.15 line spacing, including figures
(adjust to course guidelines).

The paper is the *long form* of the same story the talk and poster tell. It can carry
detail the talk can't — but it should not introduce new arguments. If something
isn't on the poster, ask whether it really belongs in the paper.

---

## Section budget

| Section | Target length | Purpose |
|---|---|---|
| Abstract | 150-200 words | Self-contained summary; the only part many readers will read. |
| Introduction | ~1 page | Motivate HARs + the open question; end with a sharp hypothesis. |
| Methods | ~1.5 pages | Reproducible. Reads like a recipe, not a story. |
| Results | ~2 pages | Three subsections: features → classification → interpretation. |
| Discussion | ~1 page | Biological interpretation, limitations, future work. |
| References | ~0.5 page | ~15 references max. |
| Figures (5) + Tables (3) | inline | See list below. |

---

## Abstract (150-200 words)

Single paragraph, structured roughly:

- **Background (1-2 sentences):** HARs are noncoding regions conserved across
  vertebrates with rapid human-lineage substitution; their distinctive functional
  contexts remain incompletely characterized.
- **Question/approach (2-3 sentences):** I compared HARs to length- and
  conservation-matched CNEs across 7 genomic features, then trained interpretable
  classifiers (logistic regression and random forest) and used SHAP to identify
  features distinguishing HARs.
- **Result (2-3 sentences):** Top-ranked feature(s); classifier AUROC; HARE5 case
  study.
- **Implication (1 sentence):** What this means for the regulatory landscape of
  human brain evolution.

Write this last.

---

## 1. Introduction (~1 page)

Three paragraphs:

**P1 — What HARs are.** Define HARs. Brief history (Pollard 2006, Capra 2013).
~96% noncoding. Multiple HAR sets in the literature; I use [the one you actually
used and why].

**P2 — Why they matter.** Enrichment near brain-expressed genes (early studies).
Specific evidence of functional consequence: Doan 2016 cognitive/behavioral
mutations, Boyd 2015 HARE5/FZD8 transgenic mouse cortex enlargement, Cui 2025
neuronal characterization. Frame these as motivating but not closing the question.

**P3 — The gap and the question.** The comparison most studies make is HARs vs
genome background, which conflates "is conserved" with "is human-accelerated". A
matched-control design isolates the HAR-specific signal. Interpretable ML lets us
ask which features carry that signal. End with the explicit hypothesis from the
proposal.

---

## 2. Methods (~1.5 pages)

Sub-headings, in order:

**2.1 HAR set.** Source, n, coordinate system, liftOver from hg19 to hg38 if used.
Cite the source paper's Table S1 or equivalent.

**2.2 Matched CNE controls.** Pre-called UCSC phastCons 100-way conserved elements
as the candidate pool. Filtering: drop any element overlapping a HAR (±10 kb buffer)
or a coding exon. Per-HAR matching on length (±25%) and mean phastCons score (±0.05).
3:1 sampling ratio. Random seed reported. QC table referenced.

**2.3 Features.** The 7 features, each with one sentence describing what it is and
how it's computed. Software versions for each external dataset (GENCODE v45, GTEx
v8, ENCODE cCRE combined track, Roadmap E081/E082 25-state imputed).

**2.4 Brain-expressed gene set.** GTEx v8 median TPM > 5 in any of 13 brain
tissues. Cite the GTEx Consortium paper.

**2.5 Classification.** Logistic regression (max_iter=5000, balanced class weight,
log1p transform on length and distance features, standardized) and random forest
(500 trees, no max depth, balanced class weight). Stratified 5-fold cross-validation.
Metrics reported: AUROC, AUPRC, accuracy, F1.

**2.6 Interpretation.** SHAP TreeExplainer on the random forest refit on all data.
Mean |SHAP| ranking and direction of effect. HARE5/2xHAR.238 case study scored
from the Boyd 2015 hg19 interval after liftOver to hg38; if absent from the
modeled Doan HAR table, it is treated as an external reference rather than a
training-set observation.

**2.7 Sensitivity analysis: brain-TSS proximity.** To test whether the reversed
brain-TSS proximity finding (HARs farther from brain-expressed TSSes than matched
CNEs) reflects a property of HARs or an artifact of conservation matching — since
conserved noncoding elements concentrate near genes — I generated a third
comparison set of random genomic intervals matched to HARs on length only
(BedTools `shuffle` with the HAR set as the source of widths), excluding HARs
and coding exons so the random set remains noncoding. Distance to the nearest
brain-expressed TSS was then compared across HAR, matched CNE, and random sets,
with pairwise Mann-Whitney U tests. Three pre-decided interpretation rules were
committed before computing the result (see `scripts/05_sensitivity_brain_tss.py`).

**2.8 Code and data availability.** GitHub repo URL, commit hash, conda env file.

---

## 3. Results (~2 pages)

Three sub-headings:

**3.1 Feature distributions, HARs vs matched CNEs.** Reference Figure 1 (the
distributions panel). Note any features showing visible univariate separation
*before* discussing the multivariate model.

**3.2 Classification performance.** Reference Figure 2 (ROC + PR curves) and
Table 1 (CV metrics by model and fold). State that performance is modest but
well above chance for both models, with consistent rankings — frame this as
support for the interpretation rather than the headline result.

**3.3 Feature importance and the HARE5 case study.** Reference Figure 3 (SHAP
beeswarm + bar) and Table 2 (top SHAP features with direction of effect). State
top-ranked feature(s) and the direction. Then introduce HARE5 (Figure 4): show
its feature values within the HAR/CNE distributions; note that it is an external
Boyd 2015 reference interval scored with the same feature-generation code.

**3.4 Sensitivity: is the brain-TSS proximity reversal a matching artifact?**
Reference Figure 5 (three-way violin of distance-to-nearest-brain-TSS) and
Table 3 (per-set median/IQR plus pairwise Mann-Whitney U). Report which of the
three pre-decided interpretations the data support: random > matched CNE > HAR
(matching artifact, strongest story); random ≈ matched CNE > HAR (the reversal
is real and not driven by matching); or random ≈ HAR < matched CNE (matched
CNEs anomalously proximal).

---

## 4. Discussion (~1 page)

**P1 — Biological interpretation.** What top SHAP features tell us about where
in the regulatory landscape human-lineage acceleration tends to fall. Emphasize
the supported claim: HARs preferentially overlap annotated regulatory elements
relative to length- and conservation-matched controls, while being farther from
gene bodies/TSSes than those controls. The sensitivity comparison against random 
length-matched intervals (§3.4) [confirms / qualifies — fill in based on the result] 
whether the distance reversal is a property of HARs or a consequence of matching 
CNEs on conservation. Connect this to long-range enhancer models and Keough 2023 
rather than claiming simple proximity to brain genes.

**P2 — Why this design works.** The matched-control + interpretable-ML combination
isolates a signal that genome-background comparisons can't, and reports it in
units (feature contributions) that are biologically interpretable.

**P3 — Limitations.** Mean phastCons remains partly definitional because the
100-way track includes human, so residual conservation signal should be framed
as expected HAR acceleration rather than an independent discovery. Annotation
bias toward studied cell types. Linear-distance proxy for regulatory targeting
(3D contacts would be better; cite Keough 2023). HAR set choice affects results
— sensitivity to alternative HAR sets is a worthwhile follow-up.

**P4 — Future work.** One paragraph. Strongest extensions: incorporate Hi-C-derived
target genes (Keough 2023), add cell-type-specific regulatory annotations from
single-cell brain atlases, replicate on the recent T2T great ape pan-genome HAR
calls (Yoo, Rhie et al. 2025).

---

## Figures

| # | Filename | What it shows | Where it's referenced |
|---|---|---|---|
| 1 | `feature_distributions.png` | Per-feature distributions, HAR vs CNE | §3.1 |
| 2 | `roc_curves.png` | ROC + PR curves, 5-fold CV mean | §3.2 |
| 3 | `shap_summary.png` | SHAP beeswarm (with `shap_bar.png` as inset) | §3.3 |
| 4 | `hare5_case_study.png` | HARE5's feature values vs HAR/CNE distributions | §3.3 |
| 5 | `sensitivity_brain_tss.png` | Distance to nearest brain-expressed TSS, 
HAR vs matched CNE vs random length-matched | §3.4 |

## Tables

| # | Filename | What it shows |
|---|---|---|
| 1 | `cv_metrics_summary.tsv` | AUROC, AUPRC, accuracy, F1 by model (mean ± SD) |
| 2 | `top_shap_features.tsv` | Ranked features by mean |SHAP|, with direction |
| 3 | `sensitivity_brain_tss.tsv` | Brain-TSS distance summary (count, mean, q25/median/q75) 
per set + pairwise Mann-Whitney U |

---

## Writing checklist

- [ ] Every method has enough detail that someone could re-run it from scratch
      using only the paper + cited public data.
- [ ] Every figure caption says what the figure *means*, not just what it shows.
- [ ] No claim in Discussion appears without a corresponding result in §3.
- [ ] Limitations section is honest, not performative — list things that actually
      could change the conclusion if addressed.
- [ ] References are uniformly formatted to course requirements.
