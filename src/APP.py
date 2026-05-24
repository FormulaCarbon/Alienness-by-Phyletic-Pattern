from pathlib import Path
import pandas as pd
import yaml
import argparse
from tqdm import tqdm
from Bio import SeqIO
import re

from genomehandler import create_raw_multifasta, create_diamond_db, run_diamond, compute_pdt, find_hgt
from dlhandler import download_url, unzip_fastas

parser = argparse.ArgumentParser()
parser.add_argument("accession")
parser.add_argument("db_path")
parser.add_argument("--config", default="config.yaml")
parser.add_argument("--blast_path", default="")
args = parser.parse_args()

with open(args.config, 'r') as configf:
    config = yaml.safe_load(configf)
    configf.close()
    
db_path = Path(args.db_path)

with open(db_path / "config.yaml", 'r') as db_conff:
    db_config = yaml.safe_load(db_conff)
    db_conff.close()

lineage_path = db_path / db_config['lineage_path']
genomes_dir = db_path / db_config['genomes_path']
blast_dir = Path(args.blast_path)

lineages = pd.read_csv(lineage_path, sep = "\t")
accession = args.accession
query = lineages[lineages["accession"] == accession].iloc[0]
url = '_'.join(query["ftp_path"].split('_')[:-1]) + "_translated_cds.faa.gz"

print("Query FASTA: " + url)
download_url(url, blast_dir / "query.faa.gz")
unzip_fastas(blast_dir)

# T1 Creation
species_T1 = lineages[lineages['species_id'] == query['species_id']]
genus_T1 = lineages[lineages["genus_id"] == query["genus_id"]]
family_T1 = lineages[lineages["family_id"] == query["family_id"]]
bucket_T1 = {
    "species": species_T1["accession"].tolist(),
    "genus": genus_T1["accession"].tolist(),
    "family": family_T1["accession"].tolist()
}

print("creating multifastas")
#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T1["species"]], blast_dir / "species_t1.fasta")
#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T1["genus"]], blast_dir / "genus_t1.fasta")
#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T1["family"]], blast_dir / "family_t1.fasta")

# T2 Creation
bucket_T2 = {}
bucket_T2["species"] = species_T1[species_T1["accession"] != accession]["accession"].tolist()
bucket_T2["genus"] = genus_T1[genus_T1["species_id"] != query['species_id']]["accession"].tolist()
bucket_T2["family"] = family_T1[family_T1["genus_id"] != query["genus_id"]]["accession"].tolist()

#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T2["species"]], blast_dir / "species_t2.fasta")
#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T2["genus"]], blast_dir / "genus_t2.fasta")
#create_raw_multifasta([genomes_dir / (acc + ".faa") for acc in bucket_T2["family"]], blast_dir / "family_t2.fasta")


# DIAMOND
# DIAMOND + PDT
groups = {
    "species_t1": bucket_T1["species"],
    "genus_t1": bucket_T1["genus"],
    "family_t1": bucket_T1["family"],
    "species_t2": bucket_T2["species"],
    "genus_t2": bucket_T2["genus"],
    "family_t2": bucket_T2["family"],
}

query_path = blast_dir / "query.faa"
query_ids = [rec.id for rec in SeqIO.parse(query_path, "fasta")]

for group, accessions in groups.items():
    level = group.split("_")[0]
    hit_path = blast_dir / "outputs" / f"{group}_hits.tsv"
    pdt_path = blast_dir / "outputs" / f"{group}_pdts.tsv"

    if len(accessions) == 0:
        # Preserve full gene list even when bucket is empty
        pd.DataFrame({
            "qseqid": query_ids,
            "genomes_with_gene": 0,
            "total_genomes": 0,
            "PDT": 0.0,
        }).to_csv(pdt_path, sep="\t", index=False)
        continue

    db_path = blast_dir / "diamond" / group
    create_diamond_db(blast_dir / f"{group}.fasta", db_path)

    run_diamond(
        query=query_path,
        db=db_path,
        out=hit_path,
        max_target_seqs=len(accessions),
        #coverage_thresh=config["coverage_threshold"],
        #identity_thresh=config["identity_thresholds"][level],
        evalue=config.get("diamond", {}).get("evalue", 10),
        threads = 16
    )

    compute_pdt(
        tsv=hit_path,
        out=pdt_path,
        query_ids=query_ids,
        expected_genomes=accessions,
        identity_thresh=config["identity_thresholds"][level],
        coverage_thresh=config["coverage_threshold"],
    )
    
df_struct = []

for group in list(groups):
    df_struct.append((str(blast_dir / "outputs" / f"{group}_pdts.tsv"), f"PDT_{group}"))

print("constructing PDT df")
# Start with the first file
df = pd.read_csv(df_struct[0][0], sep='\t', usecols=['qseqid', 'PDT'])
df = df.rename(columns={'PDT': df_struct[0][1]})

# Iteratively merge the other results
for file, colname in df_struct[1:]:
    temp = pd.read_csv(file, sep='\t', usecols=['qseqid', 'PDT'])
    temp = temp.rename(columns={'PDT': colname})
    df = pd.merge(df, temp, on='qseqid', how='outer')
df = df.fillna(0)

# Perl FLAG_pdt fallback: if species_t2 bucket is empty, classify using species_t1
if len(bucket_T2["species"]) == 0:
    df["PDT_species_for_rule"] = df["PDT_species_t1"]
else:
    df["PDT_species_for_rule"] = df["PDT_species_t2"]

types = []
modes = []
for row in tqdm(df.itertuples(), total=len(df)):
    t, m = find_hgt(row, config["pdt_cutoffs"])
    types.append(t)
    modes.append(m)

df["type"] = types
df["mode"] = modes

# Rename 'qseqid' to 'gene' if you wish
df = df.rename(columns={'qseqid': 'gene'})

print(df.head())

import re
from Bio import SeqIO

gene_coords = {}

def parse_location(location_str):
    strand = '-' if location_str.startswith('complement') else '+'
    # Remove complement( and closing ), if present
    loc = re.sub(r'^complement\((.*)\)$', r'\1', location_str)
    # Find if fuzzy edges exist
    start_fuzzy = loc.startswith('<')
    end_fuzzy = '>' in loc.split('..')[1]
    loc_nofuzzy = loc.replace('<', '').replace('>', '')
    m = re.match(r'(\d+)\.\.(\d+)', loc_nofuzzy)
    if not m:
        return None
    start, end = int(m.group(1)), int(m.group(2))
    return start, end, strand, start_fuzzy, end_fuzzy

for record in SeqIO.parse(blast_dir / "query.faa", "fasta"):
    header = record.description
    # Look for [location=...] or location=... anywhere in header
    loc_match = re.search(r'location=([^\]\s]+)', header)
    if loc_match:
        result = parse_location(loc_match.group(1))
        if result:
            start, end, strand, start_fuzzy, end_fuzzy = result
            gene_coords[record.id] = {
                "start": start,
                "end": end,
                "strand": strand,
                "start_fuzzy": start_fuzzy,
                "end_fuzzy": end_fuzzy
            }
    else:
        # Fallback: still try for any coordinates
        match = re.search(r'(\d+)\.\.(\d+)', header)
        if match:
            start, end = map(int, match.groups())
            gene_coords[record.id] = {
                "start": start,
                "end": end,
                "strand": None,
                "start_fuzzy": False,
                "end_fuzzy": False
            }

gene_pos_df = pd.DataFrame([
    {
        "gene": k,
        "start": v["start"],
        "end": v["end"],
        "strand": v["strand"],
        "start_fuzzy": v["start_fuzzy"],
        "end_fuzzy": v["end_fuzzy"]
    }
    for k, v in gene_coords.items()
])

# Merge into main table
df = df.merge(gene_pos_df, on="gene", how="left")

df.to_csv(blast_dir / "results.tsv", sep="\t", index=False)