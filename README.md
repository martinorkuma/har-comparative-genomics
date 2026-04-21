# Human Accelerated Regions (HAR) Comparatove Genomics — Project Framework

Human Accelerated Regions (HARs) in Comparative Genomics: Association of Human-Lineage Accelerated Noncoding Regions with Brain Development and Function

### Research Question

Do HARs occupy distinctive functional contexts — particularly near neurodevelopmental genes — relative to conserved non-accelerated noncoding elements, and can interpretable ML identify which features best distinguish them?

### Hypothesis

HARs will be enriched near, or linked to, genes involved in neurodevelopment and neural regulation relative to conserved non-accelerated elements and will be better distinguished by regulatory-context features than by simple sequence composition alone.

---

### Quick start

```Bash
# 1. one-time setup (creates the conda env 'har-ml' and project directories)
bash setup_env.sh

# 2. activate the env — required before every pipeline run
conda activate har-ml

# 3. Phase 1 — download HARs + reference data, liftOver, build matched CNEs
bash scripts/01_data_acquisition/run_phase1.sh
```

### Requirements

- conda or mamba
- WSL (Ubuntu) / Linux / macOS

### Pipeline phases
```Bash 
# Data acquisition 
scripts/01_data_acquisition/

# Feature engineering 
scripts/02_feature_engineering/

# Classification 
scripts/03_classification/

# Interpretation 
scripts/04_interpretation/
```

## Directory Structure

```Text
har_project/
├── setup_env.sh                       # one-command env + directory setup
├── environment/
│   ├── environment.yml                # conda dependencies (bedtools, liftOver, UCSC tools, python)
│   └── requirements.txt               # pip dependencies (pandas, sklearn, shap, ...)
├── scripts/
│   └── 01_data_acquisition/
│       ├── 00_setup_check.sh          # verifies the env before running the pipeline
│       ├── 01_download_hars.sh        # Pollard-lab HAR BEDs (hg18)
│       ├── 02_download_references.sh  # phastCons, GENCODE, liftOver chain
│       ├── 03_liftover_hars.sh        # hg18 -> hg38
│       ├── 04_build_matched_cnes.py   # length-matched CNE control set
│       ├── run_phase1.sh              # runs the whole phase
│       └── README.md                  # phase-level details
├── data/
│   ├── raw/                           # downloaded, never edited (gitignored)
│   └── processed/                     # derived by the pipeline (gitignored)
├── logs/                              # runtime logs (gitignored)
├── outputs/
│   ├── figures/                       # poster figures
│   ├── tables/                        # poster tables
│   ├── models/                        # trained classifiers
│   └── reports/                       # written paper drafts
└── .gitignore
```


## References

- Cui et al. (2025). Comparative characterization of human accelerated regions in neurons. *Nature* 640:991–999.
- Doan et al. (2016). Mutations in human accelerated regions disrupt cognition and social behavior. *Cell* 167:341–354.
- Keough et al. (2023). Three-dimensional genome rewiring in loci with human accelerated regions. *Science* 380:eabm1696.
- Levchenko et al. (2018). Human accelerated regions and other human-specific sequence variations in the context of evolution and their relevance for brain development. *Genome Biol Evol* 10:166–188.
- Pollard et al. (2006).