import argparse
from pathlib import Path
import urllib.request
from email.utils import parsedate_to_datetime
import yaml
from datetime import date
import pandas as pd
from tqdm import tqdm
import json

from refseq_preprocessor import refseq_preprocessor
from dlhandler import create_batches, download_batches, download_url
from taxahandler import check_accession_updates

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

# get db path
parser = argparse.ArgumentParser()
parser.add_argument("dbPath", help = "Path to GENOME_DB folder containing config.yaml")
parser.add_argument("-u", "--updatedb", action= "store_true")
args = parser.parse_args()

db_path = Path(args.dbPath)
conf_path = db_path / "config.yaml"

# parse config
try:
    with open(conf_path, 'r') as f:
        conf = yaml.safe_load(f)
except FileNotFoundError:
    raise FileNotFoundError("Check that the path provided has config.yaml")
    
refseq_path = db_path / conf['refseq_path']
lineage_path = db_path / conf['lineage_path']
genomes_dir = db_path / conf['genomes_path']
last_timestamp = date.fromisoformat(conf['last_updated'])
dl_cache = db_path / "temp" / "cache.json"

# get timestamp of NCBI refseq
with urllib.request.urlopen(REFSEQ_URL) as response:
    new_timestamp = parsedate_to_datetime(response.headers["Last-Modified"]).date()
    
if last_timestamp == date(1, 1, 1):
    res = input("No Database has been downloaded. Download? [y/n]: ")
    if res == 'y':
        print("Downloading RefSeq...")
        urllib.request.urlretrieve(REFSEQ_URL, refseq_path)
        conf['last_updated'] = str(new_timestamp)
        with open(conf_path, 'w') as f:
            yaml.dump(conf, f)
            
        print("Creating Lineage File...")
        refseq_preprocessor(refseq_path, lineage_path, args.updatedb)
        print("Downloading Genomes...")
        batches = create_batches(lineage_path, genomes_dir)
        dl_cache.parent.mkdir(parents = True, exist_ok = True)
        with open(dl_cache, 'w') as f:
            json.dump({"batches": batches}, f)
        download_batches(batches)
    else:
        print("Stopping Program")
        exit()
        
elif new_timestamp > last_timestamp:
    res = input("New RefSeq Found. Update Database?? [y/n]: ")
    if res == 'y':
        print("Downloading RefSeq...")
        urllib.request.urlretrieve(REFSEQ_URL, refseq_path)
        conf['last_updated'] = str(new_timestamp)
        with open(conf_path, 'w') as f:
            yaml.dump(conf, f)
        
        lineages = pd.read_csv(lineage_path, sep = '\t')
        old_accs = lineages['accession'].to_list()
        refseq = pd.read_csv(refseq_path, sep = "\t", header=0, skiprows=1, low_memory=False)
        refseq.columns = refseq.columns.str.lstrip("#").str.strip()
        new_accs = refseq["assembly_accession"].to_list()
        
        updates = check_accession_updates(old_accs, new_accs)
        
        refseq_preprocessor(refseq_path, lineage_path)
        
        # probably do not need to split this, too lazy to fix it tho
        
        if len(updates['new']) != 0:
            for update in tqdm(updates["new"], desc = "Downloading New Genomes", total = len(updates["new"])):
                org = lineages[lineages['accession'] == update]
                if len(org) > 1:
                    raise ValueError("multiple instances of same accession detected")
                
                org = org.iloc[0].to_dict()
                path = genomes_dir / org['order'] / org['family'] / org['genus'] / org['species']
                download_url(org['ftp_path'], path)
        else: 
            print("No new genomes found.")
            
        if len(updates['update']) != 0:
            for update in tqdm(updates["update"], desc = "Updating Existing Genomes", total = len(updates["update"])): 
                org = lineages[lineages['accession'] == update]
                if len(org) > 1:
                    raise ValueError("multiple instances of same accession detected")
                
                org = org.iloc[0].to_dict()
                path = genomes_dir / org['order'] / org['family'] / org['genus'] / org['species']
                download_url(org['ftp_path'], path)
        else:
            print("No genomes to update.")
        
    else:
        print("Stopping Program")
        exit()
        
else:
    print("No updates to RefSeq Found. Stopping Program")
    exit()
    



