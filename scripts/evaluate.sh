#!/usr/bin/bash

find . -type f -path '*_check/quality_report.tsv' | sed 's|^\./||; s|_check/quality_report.tsv$||' > paths.eval
python3 evaluation.py > data.eval
