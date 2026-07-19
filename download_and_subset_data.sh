#!/usr/bin/env bash
# Downloads NanoString's public CosMx SMI Lung9_Rep1 sample (2021 FFPE lung
# release, ~1.1 GB "flat data" -- counts + metadata + per-molecule
# transcript table, no images) and subsets it to a 2x2 block of fields of
# view (FOVs 1, 2, 5, 6; ~18,500 cells) so the pipeline runs in minutes
# instead of hours. See ../README.md for why Lung9_Rep1 stands in for the
# (non-public) Human Frontal Cortex dataset from the original methodology.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DATA_DIR="data"
RAW_DIR="$DATA_DIR/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images"
MINI_DIR="$DATA_DIR/mini"
TARBALL="$DATA_DIR/Lung9_Rep1_SMI_Flat_data.tar.gz"
FOVS='1|2|5|6' # 2x2 block; edit to change the FOV subset

mkdir -p "$DATA_DIR" "$MINI_DIR"

if [ ! -f "$TARBALL" ]; then
  echo "==> Downloading Lung9_Rep1 SMI Flat data (~1.1 GB)..."
  curl -sS -o "$TARBALL" \
    "https://nanostring-public-share.s3.amazonaws.com/SMI-Compressed/Lung9_Rep1/Lung9_Rep1%20SMI%20Flat%20data.tar.gz"
fi

echo "==> Extracting only the flat CSVs (skipping per-FOV TIFF/JPG images)..."
tar -xzf "$TARBALL" -C "$DATA_DIR" \
  "Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_fov_positions_file.csv" \
  "Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_metadata_file.csv" \
  "Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_exprMat_file.csv" \
  "Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_tx_file.csv"

echo "==> Subsetting to FOVs $FOVS ..."
# metadata + fov_positions use quoted fov values ("1"); exprMat/tx_file don't.
awk -F',' -v fovs="$FOVS" 'BEGIN{split(fovs,a,"|"); for(i in a) want["\""a[i]"\""]=1}
  NR==1{print; next} want[$1]{print}' \
  "$RAW_DIR/Lung9_Rep1_metadata_file.csv" > "$MINI_DIR/mini_metadata.csv"

awk -F',' -v fovs="$FOVS" 'BEGIN{split(fovs,a,"|"); for(i in a) want[a[i]]=1}
  NR==1{print; next} want[$1]{print}' \
  "$RAW_DIR/Lung9_Rep1_fov_positions_file.csv" > "$MINI_DIR/mini_fov_positions.csv"

awk -F',' -v fovs="$FOVS" 'BEGIN{split(fovs,a,"|"); for(i in a) want[a[i]]=1}
  NR==1{print; next} want[$1]{print}' \
  "$RAW_DIR/Lung9_Rep1_exprMat_file.csv" > "$MINI_DIR/mini_exprMat.csv"

awk -F',' -v fovs="$FOVS" 'BEGIN{split(fovs,a,"|"); for(i in a) want[a[i]]=1}
  NR==1{print; next} want[$1]{print}' \
  "$RAW_DIR/Lung9_Rep1_tx_file.csv" > "$MINI_DIR/mini_tx_file.csv"

echo "==> Done. Mini dataset in $MINI_DIR:"
wc -l "$MINI_DIR"/*.csv
