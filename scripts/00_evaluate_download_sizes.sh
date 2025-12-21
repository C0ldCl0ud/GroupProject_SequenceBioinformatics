#!/bin/bash

# Base data folder
DATA_DIR="../data"

# Loop through all .txt files in subfolders
find "$DATA_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.txt" | while read txtfile; do
    echo "Processing $txtfile..."

    # Get the directory of the .txt file (subfolder)
    SUBFOLDER=$(dirname "$txtfile")

    # Loop through each accession code in the file
    while read acc; do
        # Skip empty lines
        [[ -z "$acc" ]] && continue

        uid=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$acc" | xmllint --xpath 'string(//IdList/Id)' - 2>/dev/null)
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=${uid}&retmode=json" | grep -o 'total_size=\\"[0-9]*' | cut -d'"' -f2 | awk -v acc="$acc" '{s=$1; u="B KB MB GB TB"; split(u,uu); i=1; while(s>1024 && i<5){s/=1024; i++} printf "%s; %.3f; %s\n", acc, s, uu[i]}'
    done < "$txtfile"
done

acc=SRR1234567; bytes=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&term=${acc}&retmode=json" | tr -d '\n' | grep -o 'total_size="[0-9]\{1,\}"' | cut -d'"' -f2); unit="B"; val=$bytes;
