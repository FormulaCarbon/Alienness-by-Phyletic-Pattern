import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("results")
args = parser.parse_args()

results_path = args.results

results = pd.read_csv(results_path, sep="\t")

ranges = [
    (34158,87085),
    (461919,491326),
    (1932523,1961464),
    (2083238,2126304),
    (2133235,2148912),
    (1208629,1230219),
    (917453,962005),
    (868373,882872)
]

print(results.head())

out = {}

for row in results.itertuples():
    if row.type == "HGT" and 0.0 not in [row.PDT_species_t1, row.PDT_species_t2, row.PDT_genus_t1, row.PDT_genus_t2, row.PDT_family_t1, row.PDT_family_t2]:
        out[row.gene] = 0
        for r in ranges:
            if row.start >= r[0] and row.end <= r[1]:
                out[row.gene] = 1
            
vals = list(out.values())

print(sum(vals), len(vals), sum(vals)/len(vals))