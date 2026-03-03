import pandas as pd
from tqdm import tqdm

from pathlib import Path
import urllib.request
import time
import argparse

def download_url(url: str, out_path: Path, retries: int = 3, wait: float = 1.0) -> int:

    out_path.parent.mkdir(parents = True, exist_ok = True)
    
    for i in range(retries):
        try:
            with urllib.request.urlopen(url) as response:
                data = response.read()

                with open(out_path, 'wb') as file:
                    file.write(data)
            return 0
        except:
            time.sleep(wait)
    
    return 1

def download_relatives(lineages: pd.DataFrame, tax_id: int, buckets_dir: Path, strain: str = '', sample_dist: dict = {'species': 10, 'genus': 10, 'family': 10}, bucket_names: dict =  {'species': '', 'genus': '', 'family': ''},exclusive: bool = True) -> None:
    organism = lineages[(lineages['tax_id'] == tax_id) & (lineages['strain'] == strain)].iloc[0].to_dict()
    
    if exclusive:
        species_matches = lineages[(lineages['tax_id'] == tax_id) & (lineages['strain'] != strain)]
    else:
        species_matches = lineages[(lineages['tax_id'] == tax_id)]    
    species_matches = species_matches.head(sample_dist['species'])
        
    if exclusive:
        genus_matches = lineages[(lineages['genus_id'] == organism['genus_id']) & (lineages['species_id'] != organism['species_id'])]
    else:
        genus_matches = lineages[(lineages['genus_id'] == organism['genus_id'])]
    genus_matches = genus_matches.head(sample_dist['family'])
        
    if exclusive:
        family_matches = lineages[(lineages['family_id'] == organism['family_id']) & (lineages['genus_id'] != organism['genus_id'])]
    else:
        family_matches = lineages[(lineages['family_id'] == organism['family_id'])]
    family_matches = family_matches.head(sample_dist['family'])
    
    # Download Organism    
    organism_path = buckets_dir / f"TARGET_{organism['species']}_{organism['strain']}.faa.gz"
    download_url(organism['ftp_path'], organism_path)
    
    # Download Species
    species_path = buckets_dir / (f"SPECIES_{organism['species']}" if bucket_names['species'] == '' else bucket_names['species'])
    for match in tqdm(species_matches.itertuples(), desc = "Downloading Species Matches", total=sample_dist['species']):
        download_url(str(match.ftp_path), species_path / f"{str(match.genus)}_{str(match.species)}_{str(match.strain)}.faa.gz")
        
    # Download Genus
    genus_path = buckets_dir / (f"GENUS_{organism['genus']}" if bucket_names['genus'] == '' else bucket_names['genus'])
    for match in tqdm(genus_matches.itertuples(), desc = "Downloading Genus Matches", total=sample_dist['genus']):
        download_url(str(match.ftp_path), genus_path / f"{str(match.genus)}_{str(match.species)}_{str(match.strain)}.faa.gz")
        
    # Download Family
    family_path = buckets_dir / (f"FAMILY_{organism['family']}" if bucket_names['family'] == '' else bucket_names['family'])
    for match in tqdm(family_matches.itertuples(), desc = "Downloading Family Matches", total=sample_dist['family']):
        download_url(str(match.ftp_path), family_path / f"{str(match.genus)}_{str(match.species)}_{str(match.strain)}.faa.gz")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("taxId")
    parser.add_argument("lineagePath")
    parser.add_argument("bucketsPath")
    parser.add_argument("--strain", default = '', type = str)
    parser.add_argument("-e", '--exclusive', action = "store_true")
    args = parser.parse_args()
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    
    lineages = pd.read_csv(args.lineagePath, sep = "\t", header=0, low_memory=False)
    download_relatives(lineages, int(args.taxId), Path(args.bucketsPath), strain = args.strain, exclusive=args.exclusive)
