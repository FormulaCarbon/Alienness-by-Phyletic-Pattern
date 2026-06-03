import argparse
from pathlib import Path
import os
import yaml
import pandas as pd

from iohandler import header
from dlhandler import download_url, unzip_fastas
from app2helpers import *

parser = argparse.ArgumentParser()
parser.add_argument("accession", help = "RefSeq Accession (GCF_...) of query organism", type = str)
parser.add_argument("database_path", help = "Path to Genome DB", type = str)
parser.add_argument("--dir", help = "working directory", default = "/")
parser.add_argument("--config", default="config.yaml")
parser.add_argument("-m", "--mge", help="Enable Marker Gene Enrichment", action="store_true")

args = parser.parse_args()

with open(args.config, 'r') as configf:
    config = yaml.safe_load(configf)
    configf.close()
    
db_path = Path(args.database_path)

with open(db_path / "config.yaml", 'r') as db_conff:
    db_config = yaml.safe_load(db_conff)
    db_conff.close()

lineage_path = db_path / db_config['lineage_path']
genomes_dir = db_path / db_config['genomes_path']
working_dir = Path(args.dir)

lineages = pd.read_csv(lineage_path, sep = "\t")
gene_enrichment = args.mge

# Download query protein

query_acc = args.accession
query = lineages[lineages["accession"] == query_acc].iloc[0].to_dict()
print(query)
url = '_'.join(query["ftp_path"].split('_')[:-1]) + "_translated_cds.faa.gz"
header(f"Downloading {url} to {working_dir / "query.faa.gz"}...")
download_url(url, working_dir / "query.faa.gz")
unzip_fastas(working_dir)

print("Processing FASTA")
query["file_path"] = process_fasta(working_dir / "query.faa")
query["file_path"] = working_dir / "query.faa.modified"


# Header files for type 1 inclusion
acc_s1 = working_dir / "species_1_Accession_no_species.txt"
acc_g1 = working_dir / "genus_1_Accession_no_species.txt"
acc_f1 = working_dir / "family_1_Accession_no_species.txt"

# Header files for type 2 inclusion
acc_s2 = working_dir / "species_2_Accession_no_species.txt"
acc_g2 = working_dir / "genus_2_Accession_no_species.txt"
acc_f2 = working_dir / "family_2_Accession_no_species.txt"

species_id = query["species_id"]
tax_id = query["tax_id"]


flag_pdt = all_level_blast(query_acc, query, lineages, genomes_dir, working_dir)

# Analysis

# TYPE 1
input_blasthit_species = working_dir / "species_1_output_blasthits"
input_blasthit_genus = working_dir / "genus_1_output_blasthits"
input_blasthit_family = working_dir / "family_1_output_blasthits"

thresh_pident_species = config["identity_thresholds"]["species"]
thresh_pident_genus = config["identity_thresholds"]["genus"]
thresh_pident_family = config["identity_thresholds"]["family"]

if (flag_pdt == 1):
    print("processing hits for s1")
    output_blasthit_processor(
        1, 
        input_blasthit_species, 
        acc_s1, 
        thresh_pident_species,
        "species",
        query["file_path"],
        working_dir,
        config
    )
    
print("processing hits for g1")
output_blasthit_processor(
    1, 
    input_blasthit_genus, 
    acc_g1, 
    thresh_pident_genus,
    "genus",
    query["file_path"],
    working_dir,
    config
)

print("processing hits for f1")
output_blasthit_processor(
    1, 
    input_blasthit_family, 
    acc_f1, 
    thresh_pident_family,
    "family",
    query["file_path"],
    working_dir,
    config
)

# TYPE 2
input_blasthit_species = working_dir / "species_2_output_blasthits"
input_blasthit_genus = working_dir / "genus_2_output_blasthits"
input_blasthit_family = working_dir / "family_2_output_blasthits"

thresh_pident_species = config["identity_thresholds"]["species"]
thresh_pident_genus = config["identity_thresholds"]["genus"]
thresh_pident_family = config["identity_thresholds"]["family"]

if (flag_pdt == 0):
    print("processing hits for s2")
    output_blasthit_processor(
        2, 
        input_blasthit_species, 
        acc_s2, 
        thresh_pident_species,
        "species",
        query["file_path"],
        working_dir,
        config
    )
    
print("processing hits for g2")
output_blasthit_processor(
    2, 
    input_blasthit_genus, 
    acc_g2, 
    thresh_pident_genus,
    "genus",
    query["file_path"],
    working_dir,
    config
)

print("processing hits for f2")
output_blasthit_processor(
    2, 
    input_blasthit_family, 
    acc_f2, 
    thresh_pident_family,
    "family",
    query["file_path"],
    working_dir,
    config
)

calc_s1 = working_dir / f"{thresh_pident_species}_calculations_species_1.txt"
calc_g1 = working_dir / f"{thresh_pident_genus}_calculations_genus_1.txt"
calc_f1 = working_dir / f"{thresh_pident_family}_calculations_family_1.txt"

calc_s2 = working_dir / f"{thresh_pident_species}_calculations_species_2.txt"
calc_g2 = working_dir / f"{thresh_pident_genus}_calculations_genus_2.txt"
calc_f2 = working_dir / f"{thresh_pident_family}_calculations_family_2.txt"

alien_genes_finder(
    flag_pdt,
    working_dir,
    calc_s1,
    calc_g1,
    calc_f1,
    calc_s2,
    calc_g2,
    calc_f2,
    query["file_path"],
    config
)

#if (gene_enrichment):
#    marker_gene_enrichment()