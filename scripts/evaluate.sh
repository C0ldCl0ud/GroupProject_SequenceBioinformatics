#!/usr/bin/bash

find . -type d -name '*_check*' > paths.eval
python3 evaluation.py > data.eval
