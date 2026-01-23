#!/usr/bin/bash

find . -type d -name '*_check*' | sed 's|^\./||; s|_check$||' > paths.eval
python3 evaluation.py > data.eval
