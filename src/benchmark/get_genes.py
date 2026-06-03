import argparse
import gzip
import shutil
import time
import yaml
from pathlib import Path
import pandas as pd
import json
import re
import urllib.request

coord_re = re.compile(r'(\d+)\.\.(\d+)')

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

def unzip_fastas(targetPath: Path):
    for gz in targetPath.glob("*.faa.gz"): 
        out = Path(str(gz)[:-3])  # remove .gz 
        if out.exists(): 
            continue 
        print("unzipping", gz.name) 
        with gzip.open(gz, "rb") as f_in, open(out, "wb") as f_out: 
            shutil.copyfileobj(f_in, f_out) 
    

def get_span(header):
    """
    Returns (start, end) genomic span regardless of
    complement(), join(), etc.
    """
    matches = coord_re.findall(header)

    if not matches:
        return None

    coords = [(int(s), int(e)) for s, e in matches]

    start = min(s for s, e in coords)
    end = max(e for s, e in coords)

    return start, end

parser = argparse.ArgumentParser()
parser.add_argument("data", type = str)
parser.add_argument("database_path", help = "Path to Genome DB", type = str)
parser.add_argument("--dir", help = "working directory", default = "/")
parser.add_argument("--config", default="config.yaml")

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

data = {}
with open(args.data, 'r') as f:
    data = json.load(f)

alien_genes = open('alien_genes.faa', 'w')
for name in data:    
    query_acc = data[name]['accession']
    print(query_acc)
    query = lineages[lineages["accession"] == query_acc].iloc[0].to_dict()
    url = '_'.join(query["ftp_path"].split('_')[:-1]) + "_translated_cds.faa.gz"
    download_url(url, working_dir / f"{query_acc}.faa.gz")
    unzip_fastas(working_dir)
    
    with open(working_dir /  f"{query_acc}.faa") as infile:
        header = None
        seq_lines = []

        def write_if_in_range(header, seq_lines):
            if header is None:
                return

            span = get_span(header)
            if span:
                start, end = span

                for rstart, rend in data[name]['coordinates']:
                    if start <= rend and end >= rstart:  # overlap
                        fields = dict(re.findall(r'\[([^=]+)=([^\]]+)\]', header))

                        protein_id = fields.get("protein_id")
                        header = protein_id
                        alien_genes.write(f">{header}|{query_acc}\n")
                        alien_genes.writelines(seq_lines)
                        return

        for line in infile:
            if line.startswith(">"):
                write_if_in_range(header, seq_lines)
                header = line

                
                seq_lines = []
            else:
                seq_lines.append(line)
            
alien_genes.close()
    