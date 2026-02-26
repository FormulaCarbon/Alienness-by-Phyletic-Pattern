import pandas as pd

import urllib.request
import argparse
from pathlib import Path
import time
import gzip, shutil
import os

from taxahandler import search_organisms, search_organisms_by_family, get_family

def load_refseq(path: str) -> pd.DataFrame:
    refseq = pd.read_csv(path, sep = "\t", header = 0, skiprows = 1, low_memory = False)
    refseq.columns = refseq.columns.str.lstrip("#").str.strip()
    return refseq
        
def download_url(url: str, savePath: str, retries: int = 3, wait: float = 1.0) -> int:
    outPath = Path(savePath)
    outPath.parent. mkdir(parents = True, exist_ok = True)
    
    for i in range(retries):
        try:
            with urllib.request.urlopen(url) as response:
                data = response.read()

                with open(outPath, 'wb') as file:
                    file.write(data)
            return 0
        except:
            time.sleep(wait)
    
    return 1
        
def protein_fasta_url(genomeUrl: str) -> tuple[str, str]:
    url = genomeUrl.rstrip("/")
    asm = url.split('/')[-1]
    return (asm, f"{url}/{asm}_protein.faa.gz")
    
def download_relatives(refseq: pd.DataFrame, genus: str, species: str, bucketsDir: str, genus_n: int = -1, species_n: int = -1, exclusiveGenus: bool = True):
    
    genus_df = search_organisms(refseq, genus, genus_n, species if exclusiveGenus else "")
    species_df = search_organisms(refseq, species, species_n)
    
    print(f"Genus {genus}:    {len(genus_df)} items")
    print(f"Species {species}:    {len(species_df)} items")
    
    bucketsPath = Path(bucketsDir) 
    bucketsPath.mkdir(parents = True, exist_ok = True)
    
    print("Downloading genus samples")
    genusPath = bucketsPath / genus
    for row in genus_df.itertuples():
        asm, url = protein_fasta_url(row.ftp_path) # type: ignore
        print(f"Downloading {asm} | {row.organism_name}")
        download_url(url, str(genusPath / f"{asm}.faa.gz"))
    
    print("Downloading species samples")
    speciesPath = bucketsPath / species
    for row in species_df.itertuples():
        asm, url = protein_fasta_url(row.ftp_path) # type: ignore
        print(f"Downloading {asm} | {row.organism_name}")
        download_url(url, str(speciesPath / f"{asm}.faa.gz"))

def download_family(refseq: pd.DataFrame, family: str, bucketsDir: str, family_n: int = -1, excludeGenus: str = ""):
    
    family_df = search_organisms_by_family(refseq, family, family_n, excludeGenus)
    
    print(f"Family {family}:    {len(family_df)} items")
    
    bucketsPath = Path(bucketsDir) 
    bucketsPath.mkdir(parents = True, exist_ok = True)
    
    print("Downloading family samples")
    familyPath = bucketsPath / family
    for row in family_df.itertuples():
        asm, url = protein_fasta_url(row.ftp_path) # type: ignore
        print(f"Downloading {asm} | {row.organism_name}")
        download_url(url, str(familyPath / f"{asm}.faa.gz"))
    
def unzip_fastas(targetPath: Path):
    for gz in targetPath.glob("*.faa.gz"): 
        out = Path(str(gz)[:-3])  # remove .gz 
        if out.exists(): 
            continue 
        print("unzipping", gz.name) 
        with gzip.open(gz, "rb") as f_in, open(out, "wb") as f_out: 
            shutil.copyfileobj(f_in, f_out) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genus")
    parser.add_argument("species")
    parser.add_argument("refseqPath")
    parser.add_argument("bucketsPath")
    parser.add_argument("--genusSamples", default = -1, type = int)
    parser.add_argument("--speciesSamples", default = -1, type = int)
    parser.add_argument("-e", '--exclusiveGenus', action = "store_true")
    args = parser.parse_args()
    
    bucketsPath = Path(args.bucketsPath)
    shutil.rmtree(bucketsPath)
    os.mkdir(bucketsPath)
    
    refseq = load_refseq(args.refseqPath)
    pd.set_option('display.max_colwidth', None)
    download_relatives(refseq, args.genus, args.species, args.bucketsPath, genus_n = args.genusSamples, species_n = args.speciesSamples, exclusiveGenus = args.exclusiveGenus)
    
    
    unzip_fastas(bucketsPath / args.species)
    unzip_fastas(bucketsPath / args.genus)
    
    family = get_family(f"{args.genus} {args.species}")
    if family:
        download_family(refseq, family, args.bucketsPath, args.genusSamples, excludeGenus = args.genus)
        unzip_fastas(bucketsPath / family)
    else:
        print("no family found")
        
    
    
    