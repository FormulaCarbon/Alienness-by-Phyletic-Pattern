import gzip
import shutil

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

from pathlib import Path
import urllib.request
import time
import argparse
import subprocess

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
        out = Path(str(gz)[:-3])  # remove .gz 
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
                
def merge_fastas_with_genome(folder: Path, out_file: Path): # chatgpt, redo it
    with open(out_file, "w") as out:
        for f in folder.glob("*.faa"):
            genome_id = f.stem  # e.g. genus_species_strain

            for record in SeqIO.parse(f, "fasta"):
                record.id = f"{genome_id}|{record.id}"
                record.description = ""
                SeqIO.write(record, out, "fasta")
                
def make_diamond_db(genomeFile: Path, outPath: Path):
    subprocess.run([
        "diamond", "makedb",
        "--in", str(genomeFile),
        "-d", str(outPath)
    ], check=True)
    
def run_diamond(query: Path, db: Path, out: Path, threads: int = 4): # chatgpt, redo it
    subprocess.run([
        "diamond", "blastp",
        "-q", str(query),
        "-d", str(db),
        "-o", str(out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen",
        "--threads", str(threads),
        "--fast"
    ], check=True)
    
def compute_pdt(tsv: Path, identity: float, coverage: float, total_genomes: int): # chatgpt, redo it
    df = pd.read_csv(tsv, sep="\t", header=None,
        names=["qseqid","sseqid","pident","length","qlen"])

    # coverage
    df["qcov"] = df["length"] / df["qlen"]

    # filter like APP
    df = df[(df["pident"] >= identity) & (df["qcov"] >= coverage)]

    # extract genome
    df["genome"] = df["sseqid"].str.split("|").str[0]

    # compute PDT
    pdt = (
        df.groupby("qseqid")["genome"]
        .nunique()
        / total_genomes
    )

    return pdt.to_dict()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("taxId")
    parser.add_argument("lineagePath")
    parser.add_argument("bucketsPath")
    parser.add_argument("--strain", default = '', type = str)
    parser.add_argument("--familydist", default = 10, type = int)
    parser.add_argument("--genusdist", default = 10, type = int)
    parser.add_argument("--speciesdist", default = 10, type = int)
    parser.add_argument("-e", '--exclusive', action = "store_true")
    parser.add_argument("-c", '--complete', action = "store_true")
    args = parser.parse_args()
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    
    lineages = pd.read_csv(args.lineagePath, sep = "\t", header=0, low_memory=False)
    if args.complete:
        organism_path, _, _, _ =download_relatives(lineages, int(args.taxId), Path(args.bucketsPath), strain = args.strain, exclusive=False, organism_only=True)
        unzip_fastas(organism_path)
        _, t1_spec, t1_gen, t1_fam = download_relatives(lineages, int(args.taxId), Path(args.bucketsPath) / 'T1', strain = args.strain, exclusive=False, organism_skip=True, sample_dist={'species': args.speciesdist, 'genus': args.genusdist, 'family': args.familydist})  
        unzip_fastas(t1_spec) 
        unzip_fastas(t1_gen) 
        unzip_fastas(t1_fam) 
        _, t2_spec, t2_gen, t2_fam = download_relatives(lineages, int(args.taxId), Path(args.bucketsPath) / 'T2', strain = args.strain, exclusive=True, organism_skip=True, sample_dist={'species': args.speciesdist, 'genus': args.genusdist, 'family': args.familydist})
        unzip_fastas(t2_spec) 
        unzip_fastas(t2_gen) 
        unzip_fastas(t2_fam) 
    else:
        paths = download_relatives(lineages, int(args.taxId), Path(args.bucketsPath), strain = args.strain, exclusive=args.exclusive)
        for path in paths:
            unzip_fastas(path)
            
    # chatgpt, redo it
    merge_fastas_with_genome(t1_spec, Path("T1_species.faa"))
    make_diamond_db(Path("T1_species.faa"), Path("T1_species_db"))

    merge_fastas_with_genome(t1_gen, Path("T1_genus.faa"))
    make_diamond_db(Path("T1_genus.faa"), Path("T1_genus_db"))

    merge_fastas_with_genome(t1_fam, Path("T1_family.faa"))
    make_diamond_db(Path("T1_family.faa"), Path("T1_family_db"))
