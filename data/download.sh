#!/usr/bin/env bash
set -euo pipefail

INPUT_TSV=$1
DOWNLOAD_DIR=$2

if [[ ! -f "$INPUT_TSV" ]]; then
  echo "No .tsv file found."
  echo "Usage: $0 <path_to_tsv> <path_to_download_dir>"
  exit 1
fi

echo "Downloading files..."

#this goes through all available data
#needs to be adjusted to specific acession number
cut -f7 "$INPUT_TSV" | tr '\t' '\n' | while read -r url; do

    # Skip empty lines
    [[ -z "$url" ]] && continue
    cd $DOWNLOAD_DIR || 1
    echo "Downloading: $url"
    cd -
done