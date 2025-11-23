#!/usr/bin/env bash
set -euo pipefail

#This should later eventually retrieve the ftp-url for a given Accesion number
#Usage: 00_download.sh <Accession> <dataset>
#Output: echos ftp-url

FILE=$1
DATASET=$2

grep "^$FILE" "./../data/$DATASET/filereport.tsv" | cut -f7