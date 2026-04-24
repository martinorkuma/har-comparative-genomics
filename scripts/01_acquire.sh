#!/usr/bin/env bash
# 01_acquire.sh — download all data and build the HAR + matched-CNE BED files.
# ---------------------------------------------------------------------------
# Reads paths and URLs from config.yaml. Idempotent: skips downloads where the
# target file already exists.
# ---------------------------------------------------------------------------
# Outputs:
#  data/raw/*                              # raw downloads
#  data/processed/hars.hg38.bed            # HARs in hg38, BED6
#  data/processed/cnes.hg38.bed            # length-/conservation-matched CNE control set
#  data/processed/brain_genes.tsv          # gene symbol + gene_id, GTEx-brain-expressed
#  data/processed/brain_tss.hg38.bed       # TSSes of brain-expressed genes
#  data/processed/all_tss.hg38.bed         # TSSes of all GENCODE genes
#  data/processed/encode_ccre.hg38.bed     # ENCODE cCREs
#  data/processed/fetal_brain_enhancers.hg38.bed   # Roadmap E081+E082 active enhancers

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &> /dev/null && pwd)"
cd "$ROOT"

# Tiny YAML reader (only handles flat scalar keys via Python).

yget() {
    python - "$1" <<'PY'
import sys, yaml
keys = sys.argv[1].split(".")
with open("config.yaml") as fh:
    d = yaml.safe_load(fh)
for k in keys:
    if not isinstance(d, dict) or k not in d:
        sys.stderr.write(f"yget: missing key '{sys.argv[1]}'\n")
        sys.exit(1)
    d = d[k]
if isinstance(d, list):
    print("\n".join(str(x) for x in d))
else:
    print(d)
PY
}
 
LOG=logs/01_acquire.log
mkdir -p logs data/raw data/processed outputs/tables
exec > >(tee -a "$LOG") 2>&1
echo "[acquire] $(date -Iseconds) starting"
 
# Browser-ish headers so Cell Press / some CDNs don't 403 us.
UA='Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36'
 
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
 
    CHAIN_URL=$(yget hars.chain_url)
    CHAIN=data/raw/hg19ToHg38.over.chain.gz
    BED_HG19=data/raw/hars.hg19.bed
 
    download_doan() {
        local url="$1" out="$2"
        curl -L --fail --retry 3 --retry-delay 2 \
             -A "$UA" \
             -H "Accept: */*" \
             -H "Referer: https://www.cell.com/" \
             -o "$out" "$url"
    }
 
    download_pollard_fallback() {
        # Canonical 2xHAR set (Pollard/Capra), hg19 BED. Mirror on GitHub.
        # 2,649 HARs from the UCSC hgFixed mirror of Capra et al. 2013.
        local url="https://raw.githubusercontent.com/pollardlab/har_gigapipe/main/data/hars/2xHARs.hg19.bed"
        local alt="https://zoonomiaproject.org/wp-content/uploads/2021/01/2xHARs.hg19.bed"
        echo "[acquire]   fallback: fetching 2xHAR canonical set (hg19)"
        curl -L --fail --retry 3 -A "$UA" -o data/raw/hars.hg19.bed "$url" \
            || curl -L --fail --retry 3 -A "$UA" -o data/raw/hars.hg19.bed "$alt"
    }
 
    case "$HAR_SOURCE" in
        doan2016)
            URL=$(yget hars.doan2016_url)
            XLSX=data/raw/doan2016_hars.xlsx
 
            if [[ ! -s "$XLSX" ]]; then
                if ! download_doan "$URL" "$XLSX"; then
                    echo "[acquire]   WARN: Cell Press returned an error (likely 403)."
                    echo "[acquire]   Two options:"
                    echo "[acquire]     (a) download the supplement manually in a browser from"
                    echo "[acquire]         https://www.cell.com/cell/fulltext/S0092-8674(16)31169-2"
                    echo "[acquire]         and place it at $XLSX, then re-run."
                    echo "[acquire]     (b) let this script fall back to the 2xHAR canonical set now."
                    rm -f "$XLSX"
                fi
            fi
 
            if [[ -s "$XLSX" ]]; then
                python - <<'PY'
import pandas as pd
df = pd.read_excel("data/raw/doan2016_hars.xlsx")
cols = {c.lower().strip(): c for c in df.columns}
chrom = df[cols.get("chr", cols.get("chrom"))]
start = df[cols.get("start")]
end   = df[cols.get("end")]
name_col = cols.get("name", cols.get("har_id", list(df.columns)[3]))
name  = df[name_col]
chrom = chrom.astype(str).str.replace(r"^(?!chr)", "chr", regex=True)
out = pd.DataFrame({
    "chrom": chrom, "start": start.astype(int), "end": end.astype(int),
    "name": name.astype(str), "score": 0, "strand": "+",
})
out = out.dropna().sort_values(["chrom","start"])
out.to_csv("data/raw/hars.hg19.bed", sep="\t", header=False, index=False)
print(f"[acquire]   wrote {len(out)} HARs (hg19) from Doan 2016")
PY
            else
                download_pollard_fallback
            fi
            ;;
 
        pollard2xhar)
            download_pollard_fallback
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
 
    # liftOver hg19 -> hg38 (common to all source paths)
    [[ -s "$CHAIN" ]] || curl -L --fail -A "$UA" -o "$CHAIN" "$CHAIN_URL"
    UNMAPPED=data/raw/hars.unmapped.bed
    liftOver "$BED_HG19" "$CHAIN" "$HAR_HG38" "$UNMAPPED"
    n_in=$(wc -l < "$BED_HG19")
    n_out=$(wc -l < "$HAR_HG38")
    n_un=$(wc -l < "$UNMAPPED")
    echo "[acquire]   liftOver: $n_in -> $n_out mapped, $n_un unmapped"
else
    echo "[acquire] HARs already at $HAR_HG38; skipping."
fi
 
# ---------------------------------------------------------------------------
# 2. Reference genome FASTA
# ---------------------------------------------------------------------------
FASTA=data/raw/hg38.fa
if [[ ! -s "$FASTA" ]]; then
    URL=$(yget genome.fasta_url)
    echo "[acquire] downloading hg38 FASTA (~1 GB)..."
    curl -L --fail -A "$UA" -o "${FASTA}.gz" "$URL"
    gunzip "${FASTA}.gz"
fi
[[ -s "${FASTA}.fai" ]] || samtools faidx "$FASTA"
 
# ---------------------------------------------------------------------------
# 3. phastCons
# ---------------------------------------------------------------------------
PHASTCONS_BW=data/raw/hg38.phastCons100way.bw
if [[ ! -s "$PHASTCONS_BW" ]]; then
    URL=$(yget cnes.phastcons_bw_url)
    echo "[acquire] downloading phastCons 100-way bigWig (~3 GB)..."
    curl -L --fail -A "$UA" -o "$PHASTCONS_BW" "$URL"
fi
 
PHASTCONS_ELEMENTS=data/raw/hg38.phastConsElements100way.bed
if [[ ! -s "$PHASTCONS_ELEMENTS" ]]; then
    URL=$(yget cnes.phastcons_elements_url)
    RAW=data/raw/hg38.phastConsElements100way.txt.gz
    [[ -s "$RAW" ]] || curl -L --fail -A "$UA" -o "$RAW" "$URL"
    zcat "$RAW" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,"phastCons_"NR,$6}' \
                | sort -k1,1 -k2,2n > "$PHASTCONS_ELEMENTS"
    echo "[acquire]   wrote $(wc -l < "$PHASTCONS_ELEMENTS") phastCons elements"
fi
 
# ---------------------------------------------------------------------------
# 4. GENCODE + TSS
# ---------------------------------------------------------------------------
GTF=data/raw/gencode.v45.basic.annotation.gtf.gz
if [[ ! -s "$GTF" ]]; then
    URL=$(yget annotations.gencode_gtf_url)
    echo "[acquire] downloading GENCODE..."
    curl -L --fail -A "$UA" -o "$GTF" "$URL"
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
        if gtype not in ("protein_coding","lncRNA"):
            continue
        tss = start if strand == "+" else end - 1
        out.write(f"{chrom}\t{tss}\t{tss+1}\t{gid}|{gname}\t0\t{strand}\n")
out.close()
PY
    sort -k1,1 -k2,2n "$ALL_TSS" -o "$ALL_TSS"
    echo "[acquire]   wrote $(wc -l < "$ALL_TSS") TSSes"
fi
 
# ---------------------------------------------------------------------------
# 5. GTEx brain-expressed genes
# ---------------------------------------------------------------------------
BRAIN_TSS=data/processed/brain_tss.hg38.bed
if [[ ! -s "$BRAIN_TSS" ]]; then
    URL=$(yget brain_genes.gtex_median_tpm_url)
    GCT=data/raw/gtex_median_tpm.gct.gz
    [[ -s "$GCT" ]] || curl -L --fail -A "$UA" -o "$GCT" "$URL"
 
    python - <<'PY'
import yaml, pandas as pd
cfg = yaml.safe_load(open("config.yaml"))
brain_tissues = cfg["brain_genes"]["brain_tissues"]
thresh        = cfg["brain_genes"]["tpm_threshold"]
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
 
    python - <<'PY'
import pandas as pd
brain = set(pd.read_csv("data/processed/brain_genes.tsv", sep="\t")["gene_id"])
n = 0
with open("data/processed/all_tss.hg38.bed") as fh, open("data/processed/brain_tss.hg38.bed","w") as out:
    for line in fh:
        gid = line.split("\t")[3].split("|")[0]
        if gid in brain:
            out.write(line); n += 1
print(f"[acquire]   wrote {n} brain-expressed TSSes")
PY
fi
 
# ---------------------------------------------------------------------------
# 6. ENCODE cCREs
# ---------------------------------------------------------------------------
CCRE_BED=data/processed/encode_ccre.hg38.bed
if [[ ! -s "$CCRE_BED" ]]; then
    URL=$(yget annotations.encode_ccre_url)
    BB=data/raw/encodeCcreCombined.bb
    [[ -s "$BB" ]] || curl -L --fail -A "$UA" -o "$BB" "$URL"
    bigBedToBed "$BB" "$CCRE_BED"
    echo "[acquire]   wrote $(wc -l < "$CCRE_BED") cCREs"
fi
 
# ---------------------------------------------------------------------------
# 7. Roadmap fetal-brain active enhancers (hg19 -> hg38)
# ---------------------------------------------------------------------------
FBE_BED=data/processed/fetal_brain_enhancers.hg38.bed
if [[ ! -s "$FBE_BED" ]]; then
    E81_URL=$(yget annotations.roadmap_e081_url)
    E82_URL=$(yget annotations.roadmap_e082_url)
    E81=data/raw/E081_25_state.bed.gz
    E82=data/raw/E082_25_state.bed.gz
    [[ -s "$E81" ]] || curl -L --fail -A "$UA" -o "$E81" "$E81_URL"
    [[ -s "$E82" ]] || curl -L --fail -A "$UA" -o "$E82" "$E82_URL"
 
    CHAIN=data/raw/hg19ToHg38.over.chain.gz
    [[ -s "$CHAIN" ]] || curl -L --fail -A "$UA" -o "$CHAIN" "$(yget hars.chain_url)"
 
    # Read the active-enhancer state list into a temp file so Python picks it up.
    yget annotations.roadmap_active_enhancer_states > data/raw/_active_states.txt
 
    python - <<'PY'
import gzip
states = set(open("data/raw/_active_states.txt").read().splitlines())
out = open("data/raw/fetal_brain_enhancers.hg19.bed", "w")
for f in ["data/raw/E081_25_state.bed.gz", "data/raw/E082_25_state.bed.gz"]:
    sample = "E081" if "E081" in f else "E082"
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("track"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            chrom, start, end, state = parts[:4]
            if state in states:
                out.write(f"{chrom}\t{start}\t{end}\t{sample}|{state}\n")
out.close()
PY
 
    liftOver data/raw/fetal_brain_enhancers.hg19.bed "$CHAIN" \
             "$FBE_BED" data/raw/fetal_brain_enhancers.unmapped.bed
    sort -k1,1 -k2,2n "$FBE_BED" -o "$FBE_BED"
    echo "[acquire]   wrote $(wc -l < "$FBE_BED") fetal-brain enhancer intervals"
fi
 
# ---------------------------------------------------------------------------
# 8. Matched CNE control set
# ---------------------------------------------------------------------------
CNE_BED=data/processed/cnes.hg38.bed
if [[ ! -s "$CNE_BED" ]]; then
    echo "[acquire] building matched CNE control set ..."
    python scripts/_build_matched_cnes.py
fi
 
echo "[acquire] $(date -Iseconds) done."