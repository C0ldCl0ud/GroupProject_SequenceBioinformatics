#!/usr/bin/env bash
set -euo pipefail

# Aufruf: ./00_download_single.sh <dataset-name>
DATASET="${1:-}"

if [[ -z "$DATASET" ]]; then
  echo "Usage: $0 <dataset-name>"
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../data/$DATASET"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"

if [[ ! -d "$DATA_DIR" ]]; then
  echo "Fehler: Verzeichnis '$DATA_DIR' existiert nicht."
  exit 3
fi

# Funktion: prüft, ob eine Accession bereits heruntergeladen wurde.
# Annahme: prefetch legt Dateien/Ordner an, die die Accession als Präfix/Name enthalten.
# Falls prefetch anderes Muster nutzt, passe die Prüfung an.
is_downloaded() {
    local accession="$1"
    local dir="$2"
    # Suche nach Dateien/Verzeichnissen, deren Name die Accession enthält
    if find "$dir" -maxdepth 1 -name "*${accession}*" -print -quit | grep -q .; then
      return 0
    fi
    return 1
}

find "$DATA_DIR" -type f -name "*.txt" -print0 |
while IFS= read -r -d '' txtfile; do
  echo "Processing $txtfile..."
  SUBFOLDER=$(dirname "$txtfile")

  while IFS= read -r accession || [[ -n "$accession" ]]; do
  accession="${accession#"${accession%%[![:space:]]*}"}"
  accession="${accession%"${accession##*[![:space:]]}"}"
  [[ -z "$accession" ]] && continue

    if is_downloaded "$accession" "$SUBFOLDER"; then
    echo "Already downloaded: $accession (skipping)"
    continue
    fi

    echo "Prefetching $accession..."
    prefetch -O "$SUBFOLDER" -p 1 "$accession"
  done < "$txtfile"
done
