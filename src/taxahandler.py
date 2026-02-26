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
        "kingdom_id" : [],
        "phylum_id" : [],
        "class_id" : [],
        "order_id" : [],
        "family_id" : [],
        "genus_id" : [],
        "species_id": [],
        "ftp_path": ftps
    }
    
    for lineage in tqdm(lineages, desc="Creating Table"):
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(list(ranks.keys()))
        ranks = {v: k for k, v in ranks.items()}
        
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]: # slicing off the first two will exclude 'no rank' and 'domain'
            if rank in ranks:
                data[rank].append(names.get(ranks[rank], np.nan))
                data[rank + "_id"].append(ranks[rank])
            else:
                data[rank].append(np.nan)
                data[rank + "_id"].append(np.nan)
        
            
    return pd.DataFrame.from_dict(data)
        
