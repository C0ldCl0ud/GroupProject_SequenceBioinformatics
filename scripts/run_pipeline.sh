#!/usr/bin/env bash
set -euo pipefail

#orchestrates the pipeline for a given accession number

FILE=$1
DATASET=$2

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR/../data"

URL=$(bash 00_download.sh "$FILE" "$DATASET")

URL1=$(echo "$URL" | cut -d';' -f1)
URL2=$(echo "$URL" | cut -d';' -f2)

OUT_DIR="$PROJECT_DIR/$DATASET/assemblies/$FILE"
mkdir -p "$OUT_DIR"

TMP_DIR=$(mktemp -d)

#either 2 or 1 file
if [[ -n "$URL2" ]]; then
  echo "Performing short-read paired-end assembly"

  echo "Downloading $URL1 ..."
  curl -L "$URL1" | gunzip -c > "$TMP_DIR/read_1.fastq"
  echo "Downloading $URL2 ..."
  curl -L "$URL2" | gunzip -c > "$TMP_DIR/read_2.fastq"

  bash "$SCRIPT_DIR/01_assembly.sh" "$TMP_DIR/read_1.fastq" megahit "$OUT_DIR" "$TMP_DIR/read_2.fastq"

else
  echo "only one file"
fi

rm -rf "$TMP_DIR"

#curl -L "$URL" | gunzip -c | bash 01_assembly.sh $FILE