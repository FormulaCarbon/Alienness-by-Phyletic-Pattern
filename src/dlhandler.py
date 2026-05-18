from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import shutil
import os

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

from pathlib import Path
import time
import argparse
import subprocess
from itertools import batched
import re

def safe_filename(name: str) -> str:
    name = re.sub(r'[<>:"/\\\\|?*]', '_', name)
    name = name.strip('. ')
    return name[:200]

def download_url(url: str, out_path: Path, retries: int = 5):

    if out_path.exists():
        return 0

    out_path.parent.mkdir(parents=True, exist_ok=True)

    rsync_url = url.replace(
        "https://ftp.ncbi.nlm.nih.gov/",
        "rsync://ftp.ncbi.nlm.nih.gov/"
    )

    for i in range(retries):
        try:
            result = subprocess.run(
                [
                    "rsync",
                    "-a",
                    "--partial",
                    rsync_url,
                    str(out_path)
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True
            )

            if result.returncode == 0:
                return 0

            raise RuntimeError(result.stderr.strip())

        except Exception as e:
            print(f"[{i+1}/{retries}] FAILED {url}: {e}")
            time.sleep(2 ** i)

    return 1

def download_relatives(lineages: pd.DataFrame, tax_id: int, buckets_dir: Path, strain: str = '', sample_dist: dict = {'species': 10, 'genus': 10, 'family': 10}, bucket_names: dict =  {'species': '', 'genus': '', 'family': ''},exclusive: bool = True, organism_only: bool = False, organism_skip: bool = False) -> tuple:
    
    organism = lineages[(lineages['tax_id'] == tax_id) & (lineages['strain'] == strain)].iloc[0].to_dict()
    
    if exclusive:
        species_matches = lineages[(lineages['tax_id'] == tax_id) & (lineages['strain'] != strain)]
    else:
        species_matches = lineages[(lineages['tax_id'] == tax_id)]
    
    if sample_dist['species'] != 0:
        species_matches = species_matches.head(sample_dist['species'])
        
    if exclusive:
        genus_matches = lineages[(lineages['genus_id'] == organism['genus_id']) & (lineages['species_id'] != organism['species_id'])]
    else:
        genus_matches = lineages[(lineages['genus_id'] == organism['genus_id'])]
    if sample_dist['genus'] != 0:
        genus_matches = genus_matches.head(sample_dist['genus'])
        
    if exclusive:
        family_matches = lineages[(lineages['family_id'] == organism['family_id']) & (lineages['genus_id'] != organism['genus_id'])]
    else:
        family_matches = lineages[(lineages['family_id'] == organism['family_id'])]
    if sample_dist['family'] != 0:
        family_matches = family_matches.head(sample_dist['family'])
    
    # Download Organism    
    if not organism_skip:
        organism_path = buckets_dir / f"TARGET_{organism['species']}_{organism['strain']}.faa.gz"
        download_url(organism['ftp_path'], organism_path)
    
    if organism_only:
        return (organism_path, None, None, None)
    
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
    
    if organism_skip:
        return (None, species_path, genus_path, family_path)
    
    return (organism_path, species_path, genus_path, family_path)
        
def unzip_fastas(targetPath: Path):
    for gz in targetPath.glob("*.faa.gz"): 
        out = Path(str(gz)[:-3])

        if out.exists(): 
            continue 

        print("unzipping", gz.name)

        with gzip.open(gz, "rb") as f_in, open(out, "wb") as f_out: 
            shutil.copyfileobj(f_in, f_out) 
            
def merge_fastas(targetPath: Path, outPath: Path):
    with open(outPath, "wb") as f_out:
        for f in targetPath.glob("*.faa"):
            with open(f, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)
                
def merge_fastas_with_genome(folder: Path, out_file: Path):
    with open(out_file, "w") as out:
        for f in folder.glob("*.faa"):
            genome_id = f.stem

            for record in SeqIO.parse(f, "fasta"):
                record.id = f"{genome_id}|{record.id}"
                record.description = ""
                SeqIO.write(record, out, "fasta")

def create_batches(lineage_path: Path, db_path: Path, batch_size: int = 100):
    lineages = pd.read_csv(lineage_path, sep="\t", header=0, low_memory=False)

    groups = lineages.groupby('order')[[
        'accession',
        'family',
        'genus',
        'species',
        'strain',
        'ftp_path'
    ]].apply(lambda x: x.values.tolist())
    
    downloads = []
    
    for order, groups in groups.items():
        group_dir = db_path / str(order)
        
        for org in groups:
            out_path = group_dir / org[1] / org[2] / org[3]

            downloads.append((
                org[5],
                str(out_path / safe_filename(f"{org[0]}.faa.gz"))
            ))
            
    batches = list(batched(downloads, batch_size))

    return batches

def download_batches(batches, start=0, end=-1, max_workers=8):

    if end == -1:
        end = len(batches)

    print(f"Downloading from batch {start} to batch {end}")

    for i, batch in enumerate(batches[start:end]):

        print(f"Downloading batch {i} of {end-start}")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:

            futures = []

            for item in batch:
                out_path = Path(item[1])

                futures.append(
                    executor.submit(download_url, item[0], out_path)
                )

            for future in tqdm(as_completed(futures), total=len(futures)):
                future.result()

if __name__ == "__main__":
    pass