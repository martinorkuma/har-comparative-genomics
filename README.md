# Human Accelerated Regions (HAR) Comparatove Genomics — Project Framework

Human Accelerated Regions (HARs) in Comparative Genomics: Association of Human-Lineage Accelerated Noncoding Regions with Brain Development and Function

### Research Question

Do HARs occupy distinctive functional contexts — particularly near neurodevelopmental genes — relative to conserved non-accelerated noncoding elements, and can interpretable ML identify which features best distinguish them?

### Hypothesis

HARs will be enriched near, or linked to, genes involved in neurodevelopment and neural regulation relative to conserved non-accelerated elements and will be better distinguished by regulatory-context features than by simple sequence composition alone.

---

## Pipeline

### Step 1 — Data Acquisition

1. **Compile HAR set** from published datasets (Pollard lab, Keough et al. 2023, Doan et al. 2016).
2. **Construct matched comparison set** of conserved non-accelerated noncoding elements, controlling for length and conservation level.
3. **Obtain genomic/functional data** from UCSC Genome Browser, Ensembl, and NCBI Gene/GenBank.

**Tools:** UCSC Table Browser, BEDTools, NCBI datasets CLI, Ensembl BioMart.

### Step 2 — Feature Engineering

For each element (HAR or matched control), compute:

| Feature | Description |
|---|---|
| GC content | Fraction of G+C bases |
| Element length | Span in bp |
| Conservation metrics | Mean phyloP / phastCons score across the element |
| Regulatory annotation overlap | Overlap with ENCODE cCREs, enhancer marks (H3K27ac, H3K4me1) |
| Distance to nearest gene | bp to nearest TSS |
| Distance to nearest brain-related gene | bp to nearest TSS of a gene expressed in brain (BrainSpan / HPA) |
| Local TE density | Transposable element coverage in flanking window |

**Tools:** BEDTools, pandas, Biopython (SeqUtils for GC), UCSC bigWigAverageOverBed.

### Step 3 — Classification

1. **Baseline model:** Logistic regression.
2. **Nonlinear model:** Random forest.
3. **Validation:** Stratified k-fold cross-validation.
4. **Metrics:** Accuracy, precision, recall, F1, ROC-AUC.

**Tools:** scikit-learn.

### Step 4 — Interpretation

1. Feature importance from random forest.
2. SHAP values for per-feature, per-sample explanations.
3. Identify which features contribute most strongly to HAR classification.

**Tools:** SHAP, matplotlib, seaborn.


## Directory Structure

```Text
har_project/
├── data/
│   ├── hars/              # HAR BED files from published sources
│   ├── controls/          # Matched conserved non-accelerated elements
│   ├── annotations/       # cCREs, enhancer marks, TE annotations
│   ├── conservation/      # phyloP, phastCons tracks
│   └── genes/             # Gene coordinates, brain-gene lists
├── features/
│   └── feature_matrix.csv
├── models/
│   ├── logistic_reg.pkl
│   └── random_forest.pkl
├── results/
│   ├── figures/
│   ├── tables/
│   └── shap/
├── scripts/
│   ├── 01_download_data.sh
│   ├── 02_build_controls.py
│   ├── 03_compute_features.py
│   ├── 04_train_models.py
│   ├── 05_shap_analysis.py
│   └── 06_make_figures.py
├── paper/
│   └── draft.md
├── poster/
└── README.md
```


## References

- Cui et al. (2025). Comparative characterization of human accelerated regions in neurons. *Nature* 640:991–999.
- Doan et al. (2016). Mutations in human accelerated regions disrupt cognition and social behavior. *Cell* 167:341–354.
- Keough et al. (2023). Three-dimensional genome rewiring in loci with human accelerated regions. *Science* 380:eabm1696.
- Levchenko et al. (2018). Human accelerated regions and other human-specific sequence variations in the context of evolution and their relevance for brain development. *Genome Biol Evol* 10:166–188.
