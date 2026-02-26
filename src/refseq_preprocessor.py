from taxahandler import gen_taxa_db
import pandas as pd
from colorama import Fore, Back, Style, just_fix_windows_console

import argparse

just_fix_windows_console()

parser = argparse.ArgumentParser()

parser.add_argument("refseqPath")
parser.add_argument("outPath")
args = parser.parse_args()


print(Style.BRIGHT + Fore.GREEN + f"Importing RefSeq from {args.refseqPath}..." + Style.RESET_ALL)
refseq = pd.read_csv(args.refseqPath, sep = "\t", header=0, skiprows=1, low_memory=False)
refseq.columns = refseq.columns.str.lstrip("#").str.strip()
print(Style.BRIGHT + Fore.GREEN + "RefSeq Imported" + Style.RESET_ALL)
print(refseq.head())
print()


print(Style.BRIGHT + Fore.GREEN + "Filtering Bacteria Only..." + Style.RESET_ALL)
refseq = refseq[refseq["group"] == "bacteria"]
#refseq = refseq.head()
print(refseq.head())
print()

print(Style.BRIGHT + Fore.GREEN + "Adding Taxonomy Information" + Style.RESET_ALL)
taxa_db = gen_taxa_db(refseq)
print(taxa_db.head())
print()

print(Style.BRIGHT + Fore.GREEN + f"Saving to {args.outPath}..." + Style.RESET_ALL)

taxa_db.to_csv(args.outPath, sep = "\t", index = False, mode = "a")

print(Style.BRIGHT + Fore.GREEN + "Done!" + Style.RESET_ALL)
