#!/usr/bin/bash

comm -3 <(find . -type f -name 'tRNA_count.txt' -printf '%h\n' | sort -u) <(find . -type f -name 'quality_report.tsv' -path '*/_check/quality_report.tsv' -printf '%h\n' | sed 's@/_check@@' | sort -u) | sed 's|^\./||; s|_check/quality_report.tsv$||' > paths.eval
python3 evaluation.py > data.eval
