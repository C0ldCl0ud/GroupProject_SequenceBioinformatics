#!/usr/bin/env bash
set -euo pipefail

#takes in a filename(or 2 names for hybrid) and tool and performs an assembly
#Usage: 01_assembly.sh <file> <tool> <output_dir> [file2]

READ1="$1"
TOOL="$2"
OUTDIR="$3"
READ2="${4:-}"

mkdir -p "$OUTDIR"

if [[ -n "$READ2" ]]; then
  echo "Running Assembly with paired-end reads"
  $TOOL -1 "$READ1" -2 "$READ2" -t 16 -o "$OUTDIR"
else
  echo "Running Assembly with single-end reads"
  $TOOL -1 "$READ1" -t 16 -o "$OUTDIR"
fi