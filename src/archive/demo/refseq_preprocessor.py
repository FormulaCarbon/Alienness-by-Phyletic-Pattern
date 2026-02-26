from taxahandler import add_families
import pandas as pd

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("refseqPath")
parser.add_argument("outPath")
args = parser.parse_args()

print("Importing RefSeq...")
refseq = pd.read_csv(args.refseqPath, sep = "\t", header=0, skiprows=1)
refseq.columns = refseq.columns.str.lstrip("#").str.strip()

print("Filtering Bacteria Only...")
refseq = refseq[refseq["group"] == "bacteria"]

print("Adding Families...")
refseq = add_families(refseq)

with open(args.outPath, 'w') as outfile:
    outfile.write("##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.\n")

refseq.to_csv(args.outPath, sep = "\t", index = False, mode = "a")