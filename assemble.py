import sys
import subprocess
import pandas as pd
from pathlib import Path
import gzip
import shutil

if len(sys.argv) < 4:
    print("Usage: python assemble.py <datasetName> <short_reads_filename> <long_reads_filename>")
    sys.exit(1)

short_reads = "data/"+sys.argv[1]+"/"+sys.argv[2]
long_reads = "data/"+sys.argv[1]+"/"+sys.argv[3]

df_short_reads = pd.read_csv(short_reads, header=None, names=["accession"])
df_long_reads = pd.read_csv(long_reads, header=None, names=["accession"])
df_short_reads = df_short_reads.reset_index(drop=True)
df_long_reads = df_long_reads.reset_index(drop=True)

df = pd.concat([df_short_reads, df_long_reads], axis=1)
df.columns = ["short_reads", "long_reads"]

#print(df)

def dump_and_zip(col_name: str):
    #subprocess.run(
    #    ["fasterq-dump", "data/" + sys.argv[1] + "/" + row[col_name] + "/" + row[col_name] + ".sra", "-O",
    #     "data/" + sys.argv[1]], check=True)

    for fastq_file in Path("data/" + sys.argv[1]).glob("*.fastq"):
        gz_file = fastq_file.with_suffix(".fastq.gz")
        print(f"Compressing {fastq_file} -> {gz_file}")

        subprocess.run(["gzip", fastq_file], check=True)


for index, row in df.iterrows():
    print(f"Index: {index}, Accession: {row['short_reads']}, {row['long_reads']}")
    dump_and_zip('short_reads')
    dump_and_zip('long_reads')



