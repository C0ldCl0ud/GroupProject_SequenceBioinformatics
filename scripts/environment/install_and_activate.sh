#!/usr/bin/env bash

eval "$(conda shell.bash hook)"

conda env update -f environment.yml --prune
conda activate metagenomics
