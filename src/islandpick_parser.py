import pandas as pd
from ete3 import NCBITaxa

import argparse
import gzip

from taxahandler import accession2taxid

parser = argparse.ArgumentParser()
parser.add_argument("islandpick_path")

args = parser.parse_args()

islandpick = pd.read_csv(args.islandpick_path)
print('islandpick loaded')

accessions = islandpick['accession']

print('finding taxids')
taxids = [accession2taxid(a) for a in accessions]  
islandpick['accession'] = taxids

print(islandpick)


