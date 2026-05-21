from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from typing import List
from iohandler import header
from tqdm import tqdm

def create_raw_multifasta(files: List[Path], out: Path):
    out.parent.mkdir(parents = True, exist_ok = True)
    with open(out, 'w', encoding='utf-8') as outfile:
        for file in tqdm(files, total=len(files)):
            with open(file, 'r', encoding='utf-8') as infile:
                for line in infile:
                    # Sanitization - Keeps only the protein accession in the header but also adds the file the protein came from
                    if line.startswith('>'):
                        header = line[1:].split()[0]
                        outfile.write(f'>{header}|{'.'.join(file.name.split('.')[:-1])}\n')
                    else:
                        outfile.write(line)
                
def create_diamond_db(infile: Path, outfile: Path):
    outfile.parent.mkdir(parents = True, exist_ok = True)
    header(f"Creating {outfile}.dmnd from {infile}")
    subprocess.run([
        "diamond",
        "makedb",
        "--in",
        infile,
        "--db",
        outfile
    ])

def run_diamond(query: Path, db: Path, out: Path, threads: int = 4, max_target_seqs = None, coverage_thresh = None, identity_thresh = None):
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "diamond", "blastp",
        "-q", str(query),
        "-d", str(db),
        "-o", str(out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        #"--threads", str(threads),
    ]
    if max_target_seqs is not None:
        cmd.extend(["--max-target-seqs", str(max_target_seqs)])
    if coverage_thresh is not None:
        cmd.extend(["--query-cover", str(coverage_thresh)])
    if identity_thresh is not None:
        cmd.extend(["--id", str(identity_thresh)])
    subprocess.run(cmd, check=True)

def compute_pdt(tsv: Path, out: Path):
    hits = pd.read_csv(tsv, sep='\t', header=None, usecols=[0,1], names=['qseqid', 'sseqid'])
    # Extract genome accession from sseqid
    hits['genome'] = hits['sseqid'].str.split('|').str[-1]

    # Unique set of all genomes
    all_genomes = hits['genome'].unique()
    total_genomes = len(all_genomes)

    # For each qseqid, how many unique genomes does it appear in?
    pdt_df = (
        hits.groupby('qseqid')['genome']
        .nunique()
        .reset_index()
        .rename(columns={'genome': 'genomes_with_gene'})
    )
    pdt_df['total_genomes'] = total_genomes
    pdt_df['PDT'] = 100 * pdt_df['genomes_with_gene'] / total_genomes

    # Optionally write to file
    pdt_df.to_csv(out, sep='\t', index=False)
    
def find_hgt(row, pdt_thresh) -> str:
    ALIEN = "HGT"
    NATIVE = "NATIVE"
    AMB = "AMBIGUOUS"
    
    if row.PDT_species_t2 <= pdt_thresh["species"]["low"]:
        return ALIEN
    elif row.PDT_species_t2 >= pdt_thresh["species"]["high"]:
        if row.PDT_genus_t2 >= pdt_thresh["genus"]["high"]:
            if row.PDT_family_t2 >= pdt_thresh["family"]["high"]:
                return NATIVE
            elif row.PDT_family_t2 <= pdt_thresh["family"]["low"]:
                return ALIEN
            else:
                return AMB
        elif row.PDT_genus_t2 < pdt_thresh["genus"]["high"] and row.PDT_genus_t2 > pdt_thresh["species"]["low"]:
            if row.PDT_family_t1 <= pdt_thresh["family"]["low"]:
                return ALIEN
            elif row.PDT_family_t1 >= pdt_thresh["family"]["high"]:
                return NATIVE
            else:
                return AMB
        else:
            return ALIEN
    else:
        if row.PDT_family_t1 <= pdt_thresh["family"]["low"]:
            return ALIEN
        elif row.PDT_family_t1 >= pdt_thresh["family"]["high"]:
            return NATIVE
        else:
            if row.PDT_family_t1 <= pdt_thresh["family"]["low"]:
                return ALIEN
            elif row.PDT_family_t1 >= pdt_thresh["family"]["high"]:
                return NATIVE
            else:
                return AMB