# Phase 1 — Data Acquisition

This phase produces two BED files in hg38 coordinates that the rest of the
pipeline treats as the positive and negative classes for classification:

| File | Description |
|------|-------------|
| `data/processed/hars.hg38.bed` | Human Accelerated Regions (positive class) |
| `data/processed/cnes.hg38.bed` | Length-matched Conserved Non-accelerated Elements (negative class / control) |

## Scripts

Run from the project root as `bash scripts/01_data_acquisition/run_phase1.sh`,
or run each step individually in the numbered order below.

| Step | Script | What it does |
|------|--------|--------------|
| 0 | `00_setup_check.sh` | Verify `curl`, `bedtools`, `liftOver`, and Python packages are installed. Fails fast if anything is missing. |
| 1 | `01_download_hars.sh` | Download the 1,930 primate HARs and 563 mammal HARs (FDR < 0.1) from the Pollard lab website, in hg18. |
| 2 | `02_download_references.sh` | Download `phastConsElements100way` (hg38) from UCSC, GENCODE v45 basic annotation (hg38), and the hg18→hg38 liftOver chain. Derive a BED of protein-coding exons from GENCODE. |
| 3 | `03_liftover_hars.sh` | Lift both HAR sets from hg18 to hg38 using UCSC `liftOver`, then merge into a single de-duplicated union set. |
| 4 | `04_build_matched_cnes.py` | Build the CNE control set: drop phastCons elements overlapping coding exons or HARs, then for each HAR sample a length-matched candidate (±20% by default) without replacement. |

## Data sources

| Resource | Source | Assembly |
|----------|--------|----------|
| Primate HARs (1,930) | [docpollard.org — Human Acceleration in Primate Conserved Elements](https://docpollard.org/research/human-acceleration-in-primate-conserved-elements/) | hg18 |
| Mammal HARs (563) | [docpollard.org — Human Acceleration in Mammal Conserved Elements](https://docpollard.org/research/human-acceleration-in-mammal-conserved-elements/) | hg18 |
| phastConsElements100way | [hgdownload.soe.ucsc.edu](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/phastConsElements100way.txt.gz) | hg38 |
| GENCODE v45 basic | [ftp.ebi.ac.uk](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/) | hg38 |
| liftOver chain hg18→hg38 | [hgdownload.soe.ucsc.edu](https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz) | — |

## Design decisions

**Why the Pollard-lab primate set as the primary HAR source?**
It's the most widely-cited, most directly downloadable HAR set with transparent
provenance (phyloP LRT + BH FDR on 31-species multiz alignment, filtered for
synteny and against pseudogenes/repeats/segmental duplications). The smaller
mammal set is unioned in for completeness.

**Why phastConsElements100way for CNE candidates?**
It gives a large (~10⁵ elements) pool of conserved regions in the same
coordinate system as our downstream annotations. After filtering out coding
exons and HARs, whatever remains is by definition *conserved and
non-accelerated* — exactly the negative control the hypothesis calls for.

**Why match on length only at this stage?**
Pollard's HAR BED doesn't carry a phastCons LOD score directly comparable to
phastConsElements100way's LOD, so joint length+conservation matching here
would mix two different conservation scales. Phase 2 computes conservation
(mean phastCons100way and phyloP) uniformly across both sets using the bigWig
tracks, so we can verify in Phase 2 that the two groups have comparable
conservation even though we only matched on length here.

**Why sample without replacement?**
Each CNE is assigned to a single HAR and never reused, so the control set
has the same cardinality as the HAR set and there is no artificial
pseudo-replication in downstream statistics.

## Expected outputs after a clean run

```
data/
├── raw/
│   ├── 2xPrimateHARs.hg18.bed               (~1,930 rows)
│   ├── 2xMammalHARs.hg18.bed                (~563 rows)
│   ├── phastConsElements100way.txt.gz       (~10 MB)
│   ├── phastConsElements100way.hg38.bed     (~10⁵ rows)
│   ├── gencode.v45.basic.annotation.gtf.gz  (~50 MB)
│   ├── gencode.v45.coding_exons.hg38.bed    (~3–4×10⁵ merged intervals)
│   └── hg18ToHg38.over.chain.gz             (~20 MB)
└── processed/
    ├── hars.primate.hg38.bed                (most lifted, a handful drop)
    ├── hars.primate.hg38.unmapped
    ├── hars.mammal.hg38.bed
    ├── hars.mammal.hg38.unmapped
    ├── hars.hg38.bed                        (union, de-duplicated)
    └── cnes.hg38.bed                        (matched 1:1 by default)
```

Logs for each step are written to `logs/`.

## Tuning knobs

`04_build_matched_cnes.py` accepts:
- `--ratio` — CNEs per HAR (default 1). Increase to 3 or 5 if you want more
  statistical power in Phase 3 — the trade-off is class imbalance.
- `--length-tol` — fractional length tolerance for matching (default 0.2).
  Tighter values give cleaner matching but may drop HARs with unusual lengths.
- `--seed` — RNG seed for reproducibility.

## Known caveats

- Some HARs (typically <1%) fail liftOver because the hg18 region was split
  or deleted in hg38 — they appear in `*.unmapped` files. These are dropped.
- phastCons elements that are mostly-but-not-entirely inside an exon are
  dropped by `bedtools subtract -A` (conservative).
- Chromosome Y and the mitochondrial genome are excluded from both sets,
  along with alt / unplaced contigs — this keeps downstream annotation joins
  tractable.
