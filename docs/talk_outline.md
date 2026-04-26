# 10-Minute Talk Outline

**Title:** What makes a HAR a HAR? Interpretable ML on human-accelerated noncoding regions

**Total budget:** 10 minutes (~9 spoken + 1 buffer for transitions/questions during).
**Slide count target:** 8 slides. Roughly 1 minute per slide; SHAP slide gets 2.

---

## Through-line

> "I trained a small interpretable model to distinguish HARs from matched conserved
> elements, asked it which features mattered, and checked the answer at one well-known
> HAR — HARE5, which sits near *FZD8* in the developing cortex."

Every slide either sets up that sentence, delivers it, or qualifies it. If a slide
doesn't, cut it.

---

## Slide-by-slide

### Slide 1 — Title (15 s)
- Title, name, course, advisor.
- One-line subtitle: *"Which genomic features distinguish Human Accelerated Regions
  from matched conserved elements?"*

### Slide 2 — What HARs are, and the open question (90 s)
- Define HARs in one sentence: noncoding sequences conserved across vertebrates that
  show unusually rapid substitution on the human lineage.
- Why people care: ~96% noncoding, enriched near brain-expressed genes (Pollard 2006,
  Capra 2013), individual HARs implicated in cognition (Doan 2016).
- The open question: are HARs *generally* distinguished by their genomic neighborhood,
  or just sequence-level acceleration? What features matter most?
- One image: cartoon of a HAR — conserved across mouse/dog/chicken, suddenly diverged
  in human.

### Slide 3 — Approach in one diagram (60 s)
- Pipeline figure: HAR set + matched-CNE controls → 7 features → interpretable
  classifier → SHAP ranking.
- Emphasize: **"Goal is not prediction accuracy. Goal is asking the model what it learned."**
  This pre-empts the inevitable "but what's your AUC" question and frames everything
  that follows.

### Slide 4 — Matched controls + features (60 s)
- Why matched CNEs (length + conservation) instead of random genome: avoids learning
  the trivial "HARs are short conserved elements" pattern.
- The 7 features in a simple table, color-coded by category:
  - Sequence: GC content, length
  - Conservation: mean phastCons
  - Gene context: dist to nearest TSS, dist to nearest brain-expressed TSS
  - Regulatory: ENCODE cCRE overlap, fetal-brain enhancer overlap

### Slide 5 — *(THE slide.)* SHAP feature ranking (2 min)
- Big mean |SHAP| bar chart from `outputs/figures/shap_bar.png`.
- Walk through top 2-3 features. For each:
  - What it is biologically.
  - Direction of effect (HARs higher / lower than CNEs on this feature).
- This is the slide that has to be unforgettable. Take your time.
- *Don't* show the beeswarm here — it has too many moving parts for a talk audience.
  Save it for the poster where people can stare at it.

### Slide 6 — Case study: HARE5 / FZD8 (2 min)
- Motivation: Boyd 2015's HARE5/FZD8 work is the textbook example of why HARs matter.
- HARE5/2xHAR.238 is not in the modeled Doan HAR table, so I scored the lifted-over Boyd interval as an external reference using the same feature sources.
- Show `hare5_case_study.png`
- Punchline: the panel anchors the global feature ranking to the real HARE5 locus without pretending it was part of model training.

### Slide 7 — Caveats + biological interpretation (90 s)
- Three honest caveats, fast:
  1. Mean phastCons remains partly definitional because the 100-way track includes
     human; residual conservation signal should be acknowledged, not oversold.
  2. cCRE/enhancer annotations are biased toward studied cell types.
  3. Linear distance-to-gene is a weak proxy for regulatory target. I tested
     whether the surprising reversal (HARs farther from brain TSSes than matched
     CNEs) was just a matching artifact by adding a length-matched random set as
     a third group; result on the poster. 3D contacts (Keough 2023) would be a
     better proxy regardless. 
- What this means biologically: HARs are more regulatory-element-like than matched
  conserved controls, but they are not simply closer to brain-expressed TSSes.
  That makes long-range regulatory contacts the natural next hypothesis.

### Slide 8 — Summary + thanks (30 s)
- One sentence summary slide, repeating the through-line above.
- Acknowledgments: course, advisor, data sources.
- Don't ask "any questions?" — just stop. The chair will handle questions.

---

## What's *not* in the talk (deliberately)

- ROC curves and AUROC numbers → on the poster, not in the talk.
- Logistic regression results → on the poster, mentioned only as "I also trained
  a linear baseline; the ranking is consistent."
- Feature engineering details → on the poster.
- Per-feature distribution plots → on the poster.
- Sensitivity three-way comparison (HAR vs matched CNE vs random) → on the poster, 
  mentioned in one sentence on slide 7.

These all support the result; they don't *deliver* it. Putting them in the talk would
dilute the through-line.

---

## Rehearsal checklist

- [ ] Time it. Aim for 8:45 first try, you'll naturally come in around 9:30.
- [ ] Practice the SHAP slide standalone, twice. It's the hardest one to land.
- [ ] Make sure HARE5 is *named* on slide 6 and the audience knows why I picked it.
- [ ] Anticipate three questions:
  - "What's your AUC?" — *Have a number ready, but immediately reframe to interpretation.*
  - "Why these 7 features?" — *Each has a specific hypothesis (Slide 4 backup).*
  - "Why this HAR for the case study?" — *Best-validated functional HAR in the literature.*
