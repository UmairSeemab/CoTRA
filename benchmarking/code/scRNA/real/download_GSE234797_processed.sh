#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/home/umair/Downloads/CoTRA/data/scRNA/GSE234797_processed"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234797/suppl"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

files=(
  "GSE234797_CellIDs.txt.gz"
  "GSE234797_ExpressionMatrix.mtx.gz"
  "GSE234797_Genes.txt.gz"
  "GSE234797_MetaData.csv.gz"
)

for f in "${files[@]}"; do
  echo "Downloading $f"
  wget -c -O "$f" "${BASE}/${f}"
done

echo
echo "Downloaded processed GSE234797 files:"
ls -lh
