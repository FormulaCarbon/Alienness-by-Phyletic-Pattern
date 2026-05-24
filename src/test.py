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

print(results["accession"].is_unique, results['tax_id'].is_unique)
results.set_index('tax_id', verify_integrity=True)