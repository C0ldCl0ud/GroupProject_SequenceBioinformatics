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

#either 2 or 1 file
if [[ -n "$URL2" ]]; then
  echo "Performing short-read single-sample assembly"

  TMPDIR=$(mktemp -d)
  echo "Downloading $URL1 ..."
  curl -L "$URL1" | gunzip -c > "$TMPDIR/read_1.fastq"
  echo "Downloading $URL2 ..."
  curl -L "$URL2" | gunzip -c > "$TMPDIR/read_2.fastq"

  megahit -1 "$TMPDIR/read_1.fastq" -2 "$TMPDIR/read_2.fastq" -t 8 -o $PROJECT_DIR/"$DATASET"/assemblies/"$FILE"

  rm -rf "$TMPDIR"
else
  echo "only one file"
fi

#curl -L "$URL" | gunzip -c | bash 01_assembly.sh $FILE