# Deviations from the Submitted Proposal

This document records intentional changes between the submitted project
proposal (Orkuma-HAR_comparative_genomics.docx) and the implemented pipeline.
Each deviation was made to sharpen the hypothesis test, not to change it.

## 1. Feature swap: TE density → fetal-brain enhancer overlap

**Proposal (Methods, bullet 3):** "...GC content, element length, conservation
metrics, regulatory annotation overlap, distance to nearest gene, distance
to nearest brain-related gene, and local transposable element density."

**Implementation:** Local TE density is replaced with overlap of Roadmap
Epigenomics fetal-brain (E081/E082) active enhancer states.

**Rationale:** The hypothesis specifically concerns neurodevelopmental
regulatory context. Fetal-brain enhancer overlap directly probes the tissue
and developmental window most relevant to human-lineage cortical evolution
(Boyd 2015; Capra 2013). TE density primarily addresses HAR *origins* —
whether HARs derive from or sit within transposable element-rich
neighborhoods — which is a related but distinct question (Pollard 2006;
Kapusta 2013). Substituting FBE for TE density:
  - keeps the feature count at seven, as in the proposal
  - preserves the balance across feature categories (sequence,
    conservation, gene context, regulatory)
  - sharpens the regulatory category from one generic annotation (cCRE)
    to one generic + one developmentally specific annotation (FBE)

This was approved by the instructor on [DATE] before pipeline execution.

## 2. HAR source: Doan 2016 only (combines five prior HAR studies)

**Proposal:** "Compile a HAR set from published datasets..." (plural).

**Implementation:** The Doan et al. 2016 supplemental Table S1 is used as
the single HAR source. This table is itself a combined set drawn from five
prior HAR studies (Pollard 2006, Prabhakar 2006, Bird 2007, Bush & Lahn
2008, and Lindblad-Toh 2011), so the implementation effectively inherits
the multi-source design of the proposal while keeping the input
deterministic.

## 3. HARE5 case study: external reference, not training observation

**Implementation:** HARE5 (2xHAR.238) is absent from the modeled Doan HAR
table. Rather than substituting a different HAR or omitting the case study,
04_interpret.py lifts the Boyd 2015 hg19 HARE5 interval to hg38 and
computes its features from the same external sources used for HARs and
CNEs. The resulting figure places the real HARE5 within the HAR/CNE
distributions without contaminating the training set. This is documented
in scripts/04_interpret.py:build_external_hare5_row().

## 4. Added: sensitivity analysis on brain-TSS proximity

**Proposal:** No sensitivity analysis specified.

**Implementation:** `scripts/05_sensitivity_brain_tss.py` adds a three-way
comparison of distance-to-nearest-brain-expressed-TSS across HARs, matched
CNEs, and random length-matched genomic intervals. This is a follow-up to
the SHAP result that HARs sit *farther* from brain-expressed TSSes than the
matched CNEs (Spearman ≈ +0.58 between distance and SHAP), which is the
opposite of the simple "HARs near brain genes" intuition.

**Rationale:** Matched CNEs are conservation-matched, and conserved noncoding
elements concentrate near genes — so the reversal could be a property of HARs
or an artifact of the matching. Comparing both against length-only-matched
random intervals distinguishes the two. Pre-decided interpretations were
committed in the script docstring before running the analysis. The script is
not part of `run_all.sh`; it is run after the main pipeline.
