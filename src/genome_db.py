import argparse
from pathlib import Path
import yaml
from email.utils import parsedate_to_datetime
from datetime import date
import urllib.request
import os
import shutil
from dotenv import load_dotenv

import pandas as pd
from tqdm import tqdm

from refseq_preprocessor import refseq_preprocessor
from dlhandler import download_dataset
from iohandler import header, error
from dbhelpers import *

load_dotenv()

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
NCBI_API_KEY = os.getenv('NCBI_API_KEY')

parser = argparse.ArgumentParser()
parser.add_argument("dbPath")
parser.add_argument("-u", "--updatedb", action="store_true")
parser.add_argument("-m", "--max_workers", help="max workers for rehydration", default = 10)
args = parser.parse_args()

db_path = Path(args.dbPath)

conf_path = db_path / "config.yaml"

try:
    with open(conf_path, 'r') as f:
        conf = yaml.safe_load(f)
        f.close()
except FileNotFoundError as e:
    if input("No Config File Found. Create Config? [y/n]: ") == "y":
        create_config(db_path)
        with open(conf_path, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        raise e

refseq_path = db_path / conf['refseq_path']
lineage_path = db_path / conf['lineage_path']
log_path = db_path / conf["log_path"]
genomes_dir = db_path / conf['genomes_path']

last_timestamp = date.fromisoformat(conf['last_updated'])

with urllib.request.urlopen(REFSEQ_URL) as r:
    new_timestamp = parsedate_to_datetime(r.headers["Last-Modified"]).date()
    
if last_timestamp == date(1, 1, 1):
    if input("No RefSeq Found. Download Database? [y/n]: ") != "y":
        print("Stopping Program")
        exit()
        
    header("Downloading RefSeq...")
    urllib.request.urlretrieve(REFSEQ_URL, refseq_path)
    
    conf['last_updated'] = str(new_timestamp)
    with open(conf_path, 'w') as f:
        yaml.dump(conf, f)
        
    header("Creating Lineage File...")
    refseq_preprocessor(refseq_path, lineage_path, args.updatedb)
    
    Path('temp').mkdir(exist_ok=True)
    
    header("Harvesting Accessions...")
    lineages = pd.read_csv(lineage_path, sep = '\t')
    accessions = list(lineages['accession'])
    with open('temp/accessions.txt', 'w') as f:
        f.write('\n'.join(accessions))
        f.close()
    
    download_dataset(
        Path("temp/accessions.txt"), 
        Path("temp/genomes.zip"), 
        genomes_dir, 
        max_workers=int(args.max_workers),
        api_key=NCBI_API_KEY
    )
    
    fasta_dir = genomes_dir / "ncbi_dataset" / "data"
    file_count = sum(1 for x in fasta_dir.rglob('*') if x.is_file()) - 2
    with open ("log.txt", 'w') as log:
        for file in tqdm(fasta_dir.rglob("*"), total=file_count):
            try:
                if file.is_file():
                    dest_path = genomes_dir / file.name
                    if file.suffix == ".faa":
                        dest_path = genomes_dir / (file.parent.name + file.suffix)
                    shutil.move(file, dest_path)
            except Exception as e:
                log.write(str(e)+ "\n")
        log.close()
    shutil.rmtree(fasta_dir)
        
elif new_timestamp > last_timestamp:
    if input("New RefSeq Found. Update Database? [y/n]: ") != "y":
        print("Stopping Program")
        exit()
        
    new_refseq_path = db_path / ("new" + conf['refseq_path'])
    new_lineage_path =  db_path / ("new" + conf['lineage_path'])
    
    header("Downloading RefSeq...")
    urllib.request.urlretrieve(REFSEQ_URL, new_refseq_path)
    
    conf['last_updated'] = str(new_timestamp)
    with open(conf_path, 'w') as f:
        yaml.dump(conf, f)
        
    header("Creating Lineage File...")
    refseq_preprocessor(db_path / ("new" + conf['refseq_path']), new_lineage_path, args.updatedb)
    
    old_lineages = pd.read_csv(lineage_path, sep = "\t").sort_index(axis=1)
    new_lineages = pd.read_csv(new_lineage_path, sep = "\t").sort_index(axis=1)
    
    merged = old_lineages.merge(
        new_lineages,
        how="outer",
        indicator=True
    )

    removed = merged[merged["_merge"] != "left_only"]
    added = merged[merged["_merge"] != "right_only"]
    
    for row in removed.itertuples(index=False):
        os.remove(genomes_dir / row.order / row.family / row.genus / row.species / f"{row.accession}.fna" )
    
    
else:
    print("Database up to date. Stopping Program.")  
    

    
        
        