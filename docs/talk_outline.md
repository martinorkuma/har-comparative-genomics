# 10-Minute Talk Outline

**Title:** What makes a HAR a HAR? Interpretable ML on human-accelerated noncoding regions

**Total budget:** 10 minutes (~9 spoken + 1 buffer for transitions/questions during).
**Slide count target:** 8 slides. Roughly 1 minute per slide; SHAP slide gets 2.

---

## Through-line

> "We trained a small interpretable model to distinguish HARs from matched conserved
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

### Slide 6 — HARE5 / FZD8 case study (2 min)
- "What does this look like at one specific HAR we know matters?"
- HARE5 = 2xHAR.238, an enhancer near *FZD8*. Boyd et al. 2015 showed the human
  variant drives earlier and broader expression in cortical progenitors than the
  chimp variant, producing a ~12% larger neocortex in transgenic mice.
- Show `outputs/figures/hare5_case_study.png`: HARE5's value (gold marker) within the
  HAR/CNE distributions for the top 3-4 features.
- Punchline: the feature(s) the model says matter most are exactly what HARE5 exemplifies.

### Slide 7 — Caveats + biological interpretation (90 s)
- Three honest caveats, fast:
  1. Matching on length and conservation deliberately removes those signals — what's
     left is *neighborhood*, which is the question we wanted to ask.
  2. cCRE/enhancer annotations are biased toward studied cell types.
  3. Distance-to-gene is a proxy for regulatory target; 3D contacts (Keough 2023)
     would be better.
- What this means biologically: HARs aren't just sequence-accelerated dots; they sit
  in regulatory neighborhoods that look biased toward [whatever your top feature
  ends up being]. Consistent with the model that human-lineage regulatory tweaks
  near brain-relevant genes contributed to cortical evolution.

### Slide 8 — Summary + thanks (30 s)
- One sentence summary slide, repeating the through-line above.
- Acknowledgments: course, advisor, data sources.
- Don't ask "any questions?" — just stop. The chair will handle questions.

---

## What's *not* in the talk (deliberately)

- ROC curves and AUROC numbers → on the poster, not in the talk.
- Logistic regression results → on the poster, mentioned only as "we also trained
  a linear baseline; the ranking is consistent."
- Feature engineering details → on the poster.
- Per-feature distribution plots → on the poster.

These all support the result; they don't *deliver* it. Putting them in the talk would
dilute the through-line.

---

## Rehearsal checklist

- [ ] Time it. Aim for 8:45 first try, you'll naturally come in around 9:30.
- [ ] Practice the SHAP slide standalone, twice. It's the hardest one to land.
- [ ] Make sure HARE5 is *named* on slide 6 and the audience knows why we picked it.
- [ ] Anticipate three questions:
  - "What's your AUC?" — *Have a number ready, but immediately reframe to interpretation.*
  - "Why these 7 features?" — *Each has a specific hypothesis (Slide 4 backup).*
  - "Why this HAR for the case study?" — *Best-validated functional HAR in the literature.*
