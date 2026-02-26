import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data")
args = parser.parse_args()

data = pd.read_csv(args.data, sep="\t")
fracs = data[['recent', 'species_frac', 'genus_frac', 'family_frac']]