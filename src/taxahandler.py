import argparse

from pathlib import Path
import pandas as pd
from ete3 import NCBITaxa
from tqdm import tqdm
from colorama import Fore, Back, Style
import numpy as np

def gen_taxa_db(target: pd.DataFrame) -> pd.DataFrame:
    ncbi = NCBITaxa()
    
    print("Cleaning Data")
    clean = target.dropna(subset = ["species_taxid", "ftp_path"])
    
    taxIDs = clean["species_taxid"].astype(int).tolist()
    ftps = clean['ftp_path'].astype(str).tolist()
    
    # convert to fasta link
    urls = [f"{genomeUrl.rstrip("/")}/{genomeUrl.rstrip("/").split('/')[-1]}_protein.faa.gz" for genomeUrl in ftps]
    
    strain = [items.split(', /')[0].split('=')[1] if items != 'na' else '' for items in clean['infraspecific_name'].astype(str).tolist()]
    lineages = [ncbi.get_lineage(taxID) for taxID in tqdm(taxIDs, desc="Getting Lineage Information")]
    
    data = {
        "tax_id": taxIDs,
        "kingdom" : [],
        "phylum" : [],
        "class" : [],
        "order" : [],
        "family" : [],
        "genus" : [],
        "species": [],
        "strain": strain,
        "kingdom_id" : [],
        "phylum_id" : [],
        "class_id" : [],
        "order_id" : [],
        "family_id" : [],
        "genus_id" : [],
        "species_id": [],
        "ftp_path": urls
    }
    
    for lineage in tqdm(lineages, desc="Creating Table"):
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(list(ranks.keys()))
        ranks = {v: k for k, v in ranks.items()}
        
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]: # slicing off the first two will exclude 'no rank' and 'domain'
            if rank in ranks:
                data[rank].append(names.get(ranks[rank], np.nan))
                data[rank + "_id"].append(int(ranks[rank]))
            else:
                data[rank].append(np.nan)
                data[rank + "_id"].append(np.nan)
        
    return pd.DataFrame.from_dict(data)

def graph(lineages: pd.DataFrame, outfile: Path):
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(lineages['tax_id'].astype(str).tolist())
    
    for node in tree.traverse():
        if node.name.isdigit():
            names = ncbi.get_taxid_translator([int(node.name)])
            if names:
                node.name = names[int(node.name)]
    
    tree.write(outfile=str(outfile))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("lineagePath")
    parser.add_argument("outFile")
    args = parser.parse_args()
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)

    
    lineages = pd.read_csv(args.lineagePath, sep = "\t", header=0, low_memory=False)
    graph(lineages, Path(args.outFile))
        
