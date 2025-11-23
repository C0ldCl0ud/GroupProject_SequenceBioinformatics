import sys
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: python assemble.py <path_to_short_reads> <path_to_long_reads>")
    sys.exit(1)

short_reads = sys.argv[1]
long_reads = sys.argv[2]

df_short_reads = pd.read_csv(short_reads, header=None, names=["accession"])
df_long_reads = pd.read_csv(long_reads, header=None, names=["accession"])
df_short_reads = df_short_reads.reset_index(drop=True)
df_long_reads = df_long_reads.reset_index(drop=True)

df = pd.concat([df_short_reads, df_long_reads], axis=1)
df.columns = ["short_reads", "long_reads"]

#print(df)

for index, row in df.iterrows():
    print(f"Index: {index}, Accession: {row['short_reads']}, {row['long_reads']}")




