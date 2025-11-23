#!/usr/bin/env bash

# Base data folder
DATA_DIR="data"

# Loop through all .txt files in subfolders
find "$DATA_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.txt" | while read txtfile; do
    echo "Processing $txtfile..."

    # Get the directory of the .txt file (subfolder)
    SUBFOLDER=$(dirname "$txtfile")

    # Loop through each accession code in the file
    while read accession; do
        # Skip empty lines
        [[ -z "$accession" ]] && continue

        echo "Prefetching $accession..."
        prefetch -O "$SUBFOLDER" "$accession"
    done < "$txtfile"
done