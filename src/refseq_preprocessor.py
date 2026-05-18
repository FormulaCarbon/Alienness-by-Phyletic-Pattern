from taxahandler import gen_taxa_db
import pandas as pd
from colorama import Fore, Back, Style, just_fix_windows_console

import argparse


def refseq_preprocessor(refseqPath, outPath, update = False):
    just_fix_windows_console()

    print(Style.BRIGHT + Fore.GREEN + f"Importing RefSeq from {refseqPath}..." + Style.RESET_ALL)
    refseq = pd.read_csv(refseqPath, sep = "\t", header=0, skiprows=1, low_memory=False)
    refseq.columns = refseq.columns.str.lstrip("#").str.strip()
    print(Style.BRIGHT + Fore.GREEN + "RefSeq Imported" + Style.RESET_ALL)
    print(refseq.head())
    print()


    print(Style.BRIGHT + Fore.GREEN + "Filtering Bacteria Only..." + Style.RESET_ALL)
    refseq = refseq[(refseq["group"] == "bacteria") & (refseq["version_status"] == "latest") & (refseq["assembly_level"] == "Complete Genome")]
    refseq = refseq.loc[refseq.groupby(['assembly_accession', 'infraspecific_name'])['seq_rel_date'].idxmax()]
    #refseq = refseq.head()
    print(refseq.head())
    print()

    print(Style.BRIGHT + Fore.GREEN + "Adding Taxonomy Information" + Style.RESET_ALL)
    taxa_db = gen_taxa_db(refseq, update)
    print(taxa_db.head())
    print()

    print(Style.BRIGHT + Fore.GREEN + f"Saving to {outPath}..." + Style.RESET_ALL)

    taxa_db.to_csv(outPath, sep = "\t", index = False, mode = "w", float_format='%.0f')

    print(Style.BRIGHT + Fore.GREEN + "Done!" + Style.RESET_ALL)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("refseq_path")
    parser.add_argument("out_path")

    args = parser.parse_args()
    refseq_preprocessor(args.refseq_path, args.out_path)
