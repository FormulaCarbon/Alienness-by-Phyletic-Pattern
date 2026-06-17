from taxahandler import gen_taxa_db
import pandas as pd
from colorama import Fore, Back, Style, just_fix_windows_console
import numpy as np

import argparse

just_fix_windows_console()

parser = argparse.ArgumentParser()

parser.add_argument("lineagePath")
parser.add_argument("outPath")
args = parser.parse_args()


def split_strain(i):
    try:
        name = i.split('_strain-')[0]
        strain = i.split('_strain-')[1]
    except IndexError:
        name = name
        strain = name
    return (name, strain)

print(Style.BRIGHT + Fore.GREEN + f"Importing Soham's Lineage File from {args.lineagePath}..." + Style.RESET_ALL)
lineage = pd.read_csv(args.lineagePath, sep = "\t", header=0, low_memory=False)
print(Style.BRIGHT + Fore.GREEN + "Imported" + Style.RESET_ALL)
print(lineage.head())
print()


print(Style.BRIGHT + Fore.GREEN + "Filtering Bacteria Only..." + Style.RESET_ALL)
lineage = lineage[lineage["kingdom"] == "Bacteria"]
#lineage = lineage.head()
print(lineage.head())
print()

print(Style.BRIGHT + Fore.GREEN + "Reformatting" + Style.RESET_ALL)
accs = ['_'.join(genomeUrl.split('/')[-1].split('_')[:2]) for genomeUrl in lineage['FTP_link'].astype(str).tolist()]
urls = [f"https:{genomeUrl.split(":")[1]}/{genomeUrl.rstrip("/").split('/')[-1]}_protein.faa.gz" for genomeUrl in lineage['FTP_link'].astype(str).tolist()]
data = {
        "accession": accs,
        "tax_id": lineage['TAX_ID'].astype(int).tolist(),
        "kingdom" : lineage['kingdom'].astype(str).tolist(),
        "phylum" : lineage['phylum'].astype(str).tolist(),
        "class" : lineage['class'].astype(str).tolist(),
        "order" : lineage['order'].astype(str).tolist(),
        "family" : lineage['family'].astype(str).tolist(),
        "genus" : lineage['genus'].astype(str).tolist(),
        "species": [split_strain(i)[0] for i in lineage['Species_name'].astype(str).tolist()],
        "strain": [split_strain(i)[1] for i in lineage['Species_name'].astype(str).tolist()],
        "kingdom_id" : lineage['kingdom.1'].astype(int).tolist(),
        "phylum_id" : lineage['phylum.1'].astype(int).tolist(),
        "class_id" : lineage['class.1'].astype(int).tolist(),
        "order_id" : lineage['order.1'].astype(int).tolist(),
        "family_id" : lineage['family.1'].astype(int).tolist(),
        "genus_id" : lineage['genus.1'].astype(int).tolist(),
        "species_id": lineage['SPECIES_ID'].astype(int).tolist(),
        "ftp_path": urls
    }
print(urls[:5])
df = pd.DataFrame.from_dict(data)
lineages = pd.read_csv("H:/GENOME_DB/BACTERIAL_LINEAGES.txt", sep="\t")
taxa_db = lineages[lineages["accession"].isin(df["accession"])]
print(taxa_db.head())
print()

print(Style.BRIGHT + Fore.GREEN + f"Saving to {args.outPath}..." + Style.RESET_ALL)

df.to_csv(args.outPath, sep = "\t", index = False, mode = "w", float_format='%.0f')

print(Style.BRIGHT + Fore.GREEN + "Done!" + Style.RESET_ALL)
