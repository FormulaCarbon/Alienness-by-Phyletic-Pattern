import pandas as pd
from ete3 import NCBITaxa
from tqdm import tqdm

import argparse
import gzip
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("islandpick_path")
parser.add_argument("map_path")
parser.add_argument("bacterial_lineages_path")
args = parser.parse_args()

islandpick = pd.read_csv(args.islandpick_path)
print('islandpick loaded')

mappings = pd.read_csv(args.map_path, sep = "\t").set_index('accession')['taxid'].to_dict()

lineages = pd.read_csv(args.bacterial_lineages_path, sep = "\t")

accessions = islandpick['accession']
counts = Counter(accessions)

taxids = []
remove = []
for i, acc in enumerate(islandpick['accession']):
    try:
        taxids.append(mappings[acc[:-2]])
    except KeyError:
        remove.append(i)

sort_taxids = list(Counter(taxids).keys())
sort_taxids.sort()
        
print(sort_taxids[:10])
        
filtered = lineages[lineages['tax_id'].isin(taxids)]
filtered.to_csv('islandpick_lineages.txt', sep = '\t', index = False)
        

islandpick = islandpick.drop(remove)

print(len(islandpick), len(remove), len(taxids))

islandpick['accession'] = taxids
islandpick = islandpick.rename({'accession': 'taxid'})
islandpick.to_csv('islandpick_taxids.txt', sep = '\t',  index = False)




