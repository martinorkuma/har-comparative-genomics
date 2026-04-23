#!/usr/bin/env bash
# 01_acquire.sh — download all data and build the HAR + matched-CNE BED files.
#
# Reads paths and URLs from config.yaml. Idempotent: skips downloads where the
# target file already exists.
#
# Outputs:
#   data/raw/*                              # raw downloads
#   data/processed/hars.hg38.bed            # HARs in hg38, BED6
#   data/processed/cnes.hg38.bed            # length-/conservation-matched CNE control set
#   data/processed/brain_genes.tsv          # gene symbol + gene_id, GTEx-brain-expressed
#   data/processed/brain_tss.hg38.bed       # TSSes of brain-expressed genes
#   data/processed/all_tss.hg38.bed         # TSSes of all GENCODE genes
#   data/processed/encode_ccre.hg38.bed     # ENCODE cCREs
#   data/processed/fetal_brain_enhancers.hg38.bed   # Roadmap E081+E082 active enhancers

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &> /dev/null && pwd)"
cd "$ROOT"

# Tiny YAML reader (only handles flat scalar keys via Python).
yget() {
    python -c "import yaml,sys; d=yaml.safe_load(open('config.yaml')); \
keys='$1'.split('.'); v=d; \
[v := v[k] for k in keys]; \
print(v if not isinstance(v, list) else '\n'.join(v))"
}

LOG=logs/01_acquire.log
exec > >(tee -a "$LOG") 2>&1
echo "[acquire] $(date -Iseconds) starting"

mkdir -p data/raw data/processed

# ---------------------------------------------------------------------------
# 1. HARs
# ---------------------------------------------------------------------------
HAR_HG38=data/processed/hars.hg38.bed

if [[ -s data/raw/hars.user.bed ]]; then
    echo "[acquire] using user-supplied HAR BED at data/raw/hars.user.bed"
    cp data/raw/hars.user.bed "$HAR_HG38"
elif [[ ! -s "$HAR_HG38" ]]; then
    HAR_SOURCE=$(yget hars.source)
    echo "[acquire] HAR source = $HAR_SOURCE"
    case "$HAR_SOURCE" in
        doan2016)
            URL=$(yget hars.doan2016_url)
            XLSX=data/raw/doan2016_hars.xlsx
            BED_HG19=data/raw/hars.hg19.bed
            [[ -s "$XLSX" ]] || curl -L --fail -o "$XLSX" "$URL"

            # Convert the Doan supplemental table to BED6.
            # NOTE: column names vary by paper supplement; verify with `python -c "import pandas as pd; print(pd.read_excel('data/raw/doan2016_hars.xlsx').head())"`
            python - <<'PY'
import pandas as pd
df = pd.read_excel("data/raw/doan2016_hars.xlsx")
# Standardize column names — the Doan supplement uses 'Chr', 'Start', 'End', 'Name'.
# Adjust here if your downloaded version differs.
cols = {c.lower().strip(): c for c in df.columns}
chrom = df[cols.get("chr", cols.get("chrom"))]
start = df[cols.get("start")]
end   = df[cols.get("end")]
name  = df[cols.get("name", cols.get("har_id", list(df.columns)[3]))]
# enforce 'chr' prefix
chrom = chrom.astype(str).str.replace(r"^(?!chr)", "chr", regex=True)
out = pd.DataFrame({
    "chrom": chrom, "start": start.astype(int), "end": end.astype(int),
    "name": name.astype(str), "score": 0, "strand": "+",
})
out = out.dropna().sort_values(["chrom","start"])
out.to_csv("data/raw/hars.hg19.bed", sep="\t", header=False, index=False)
print(f"[acquire]   wrote {len(out)} HARs (hg19) -> data/raw/hars.hg19.bed")
PY

            # liftOver hg19 -> hg38
            CHAIN_URL=$(yget hars.chain_url)
            CHAIN=data/raw/hg19ToHg38.over.chain.gz
            [[ -s "$CHAIN" ]] || curl -L --fail -o "$CHAIN" "$CHAIN_URL"
            UNMAPPED=data/raw/hars.unmapped.bed
            liftOver "$BED_HG19" "$CHAIN" "$HAR_HG38" "$UNMAPPED"
            n_in=$(wc -l < "$BED_HG19")
            n_out=$(wc -l < "$HAR_HG38")
            n_un=$(wc -l < "$UNMAPPED")
            echo "[acquire]   liftOver: $n_in -> $n_out mapped, $n_un unmapped"
            ;;
        capra2013)
            echo "ERROR: capra2013 source not yet wired; supply via data/raw/hars.user.bed instead." >&2
            exit 1
            ;;
        user)
            echo "ERROR: hars.source=user but data/raw/hars.user.bed missing." >&2
            exit 1
            ;;
        *)
            echo "ERROR: unknown hars.source '$HAR_SOURCE'." >&2
            exit 1
            ;;
    esac
else
    echo "[acquire] HARs already at $HAR_HG38; skipping."
fi

# ---------------------------------------------------------------------------
# 2. Reference genome FASTA (only needed for GC content)
# ---------------------------------------------------------------------------
FASTA=data/raw/hg38.fa
if [[ ! -s "$FASTA" ]]; then
    URL=$(yget genome.fasta_url)
    echo "[acquire] downloading hg38 FASTA (~1 GB)..."
    curl -L --fail -o "${FASTA}.gz" "$URL"
    gunzip "${FASTA}.gz"
fi
[[ -s "${FASTA}.fai" ]] || samtools faidx "$FASTA"

# ---------------------------------------------------------------------------
# 3. phastCons — both pre-called elements (CNE candidate pool) AND the bigwig
#    (per-base scores, used as a feature on each HAR/CNE).
# ---------------------------------------------------------------------------
PHASTCONS_BW=data/raw/hg38.phastCons100way.bw
if [[ ! -s "$PHASTCONS_BW" ]]; then
    URL=$(yget cnes.phastcons_bw_url)
    echo "[acquire] downloading phastCons 100-way bigWig (~3 GB)..."
    curl -L --fail -o "$PHASTCONS_BW" "$URL"
fi

PHASTCONS_ELEMENTS=data/raw/hg38.phastConsElements100way.bed
if [[ ! -s "$PHASTCONS_ELEMENTS" ]]; then
    URL=$(yget cnes.phastcons_elements_url)
    RAW=data/raw/hg38.phastConsElements100way.txt.gz
    [[ -s "$RAW" ]] || curl -L --fail -o "$RAW" "$URL"
    # The UCSC text dump has columns: bin, chrom, chromStart, chromEnd, lod, score
    # Strip 'bin' and emit BED5.
    zcat "$RAW" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,"phastCons_"NR,$6}' \
                | sort -k1,1 -k2,2n > "$PHASTCONS_ELEMENTS"
    echo "[acquire]   wrote $(wc -l < "$PHASTCONS_ELEMENTS") phastCons elements -> $PHASTCONS_ELEMENTS"
fi

# ---------------------------------------------------------------------------
# 4. GENCODE annotation + TSS BEDs
# ---------------------------------------------------------------------------
GTF=data/raw/gencode.v45.basic.annotation.gtf.gz
if [[ ! -s "$GTF" ]]; then
    URL=$(yget annotations.gencode_gtf_url)
    echo "[acquire] downloading GENCODE..."
    curl -L --fail -o "$GTF" "$URL"
fi

ALL_TSS=data/processed/all_tss.hg38.bed
if [[ ! -s "$ALL_TSS" ]]; then
    echo "[acquire] extracting TSSes from GENCODE..."
    python - <<'PY'
import gzip, re
out = open("data/processed/all_tss.hg38.bed", "w")
with gzip.open("data/raw/gencode.v45.basic.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene": continue
        chrom, start, end, strand = f[0], int(f[3])-1, int(f[4]), f[6]
        attrs = dict(re.findall(r'(\w+) "([^"]+)"', f[8]))
        gid = attrs.get("gene_id","").split(".")[0]
        gname = attrs.get("gene_name","")
        gtype = attrs.get("gene_type","")
        if gtype not in ("protein_coding","lncRNA"):  # keep it focused
            continue
        tss = start if strand == "+" else end - 1
        out.write(f"{chrom}\t{tss}\t{tss+1}\t{gid}|{gname}\t0\t{strand}\n")
out.close()
PY
    sort -k1,1 -k2,2n "$ALL_TSS" -o "$ALL_TSS"
    echo "[acquire]   wrote $(wc -l < "$ALL_TSS") TSSes -> $ALL_TSS"
fi

# ---------------------------------------------------------------------------
# 5. GTEx brain-expressed gene list -> brain_tss.hg38.bed
# ---------------------------------------------------------------------------
BRAIN_GENES=data/processed/brain_genes.tsv
BRAIN_TSS=data/processed/brain_tss.hg38.bed
if [[ ! -s "$BRAIN_TSS" ]]; then
    URL=$(yget brain_genes.gtex_median_tpm_url)
    GCT=data/raw/gtex_median_tpm.gct.gz
    [[ -s "$GCT" ]] || curl -L --fail -o "$GCT" "$URL"

    python - <<'PY'
import gzip, yaml, pandas as pd
cfg = yaml.safe_load(open("config.yaml"))
brain_tissues = cfg["brain_genes"]["brain_tissues"]
thresh        = cfg["brain_genes"]["tpm_threshold"]

# GTEx GCT: 2 header lines, then a table whose first two cols are gene_id, Description.
df = pd.read_csv("data/raw/gtex_median_tpm.gct.gz", sep="\t", skiprows=2, compression="gzip")
df["gene_id_short"] = df["Name"].str.split(".").str[0]
have = [t for t in brain_tissues if t in df.columns]
missing = set(brain_tissues) - set(have)
if missing:
    print("[acquire]   warning: GTEx columns not found:", missing)
brain_max = df[have].max(axis=1)
brain_set = df.loc[brain_max > thresh, ["gene_id_short", "Description"]]
brain_set.columns = ["gene_id", "gene_name"]
brain_set.to_csv("data/processed/brain_genes.tsv", sep="\t", index=False)
print(f"[acquire]   {len(brain_set)} genes with brain median TPM > {thresh}")
PY

    # Filter all_tss.hg38.bed to brain genes
    python - <<'PY'
import pandas as pd
brain = set(pd.read_csv("data/processed/brain_genes.tsv", sep="\t")["gene_id"])
with open("data/processed/all_tss.hg38.bed") as fh, open("data/processed/brain_tss.hg38.bed","w") as out:
    n = 0
    for line in fh:
        gid = line.split("\t")[3].split("|")[0]
        if gid in brain:
            out.write(line); n += 1
print(f"[acquire]   wrote {n} brain-expressed TSSes -> data/processed/brain_tss.hg38.bed")
PY
fi

# ---------------------------------------------------------------------------
# 6. ENCODE cCREs (any-cell-type combined)
# ---------------------------------------------------------------------------
CCRE_BED=data/processed/encode_ccre.hg38.bed
if [[ ! -s "$CCRE_BED" ]]; then
    URL=$(yget annotations.encode_ccre_url)
    BB=data/raw/encodeCcreCombined.bb
    [[ -s "$BB" ]] || curl -L --fail -o "$BB" "$URL"
    bigBedToBed "$BB" "$CCRE_BED"
    echo "[acquire]   wrote $(wc -l < "$CCRE_BED") cCREs -> $CCRE_BED"
fi

# ---------------------------------------------------------------------------
# 7. Roadmap fetal-brain active enhancers (E081 + E082, 25-state model)
# ---------------------------------------------------------------------------
FBE_BED=data/processed/fetal_brain_enhancers.hg38.bed
if [[ ! -s "$FBE_BED" ]]; then
    E81_URL=$(yget annotations.roadmap_e081_url)
    E82_URL=$(yget annotations.roadmap_e082_url)
    E81=data/raw/E081_25_state.bed.gz
    E82=data/raw/E082_25_state.bed.gz
    [[ -s "$E81" ]] || curl -L --fail -o "$E81" "$E81_URL"
    [[ -s "$E82" ]] || curl -L --fail -o "$E82" "$E82_URL"

    # Roadmap chromHMM is in hg19 — needs liftOver to hg38.
    CHAIN=data/raw/hg19ToHg38.over.chain.gz
    [[ -s "$CHAIN" ]] || curl -L --fail -o "$CHAIN" "$(yget hars.chain_url)"

    # Active enhancer states from config
    STATES_RAW=$(yget annotations.roadmap_active_enhancer_states)
    STATES_REGEX=$(echo "$STATES_RAW" | paste -sd'|' -)

    python - <<PY
import gzip, re
states = set("$STATES_RAW".split("\n"))
out = open("data/raw/fetal_brain_enhancers.hg19.bed", "w")
for f in ["data/raw/E081_25_state.bed.gz", "data/raw/E082_25_state.bed.gz"]:
    sample = "E081" if "E081" in f else "E082"
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("track"): continue
            chrom, start, end, state = line.rstrip("\n").split("\t")[:4]
            if state in states:
                out.write(f"{chrom}\t{start}\t{end}\t{sample}|{state}\n")
out.close()
PY

    liftOver data/raw/fetal_brain_enhancers.hg19.bed "$CHAIN" \
             "$FBE_BED" data/raw/fetal_brain_enhancers.unmapped.bed
    sort -k1,1 -k2,2n "$FBE_BED" -o "$FBE_BED"
    echo "[acquire]   wrote $(wc -l < "$FBE_BED") fetal-brain enhancer intervals -> $FBE_BED"
fi

# ---------------------------------------------------------------------------
# 8. Build matched CNE control set
# ---------------------------------------------------------------------------
CNE_BED=data/processed/cnes.hg38.bed
if [[ ! -s "$CNE_BED" ]]; then
    echo "[acquire] building matched CNE control set ..."
    python scripts/_build_matched_cnes.py
fi

echo "[acquire] $(date -Iseconds) done."
