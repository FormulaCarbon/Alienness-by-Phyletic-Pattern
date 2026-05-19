from taxahandler import gen_taxa_db
import pandas as pd
from iohandler import header

import argparse


def refseq_preprocessor(refseqPath, outPath, update = False):

    header(f"Importing RefSeq from {refseqPath}...")
    refseq = pd.read_csv(refseqPath, sep = "\t", header=0, skiprows=1, low_memory=False)
    refseq.columns = refseq.columns.str.lstrip("#").str.strip()
    header("RefSeq Imported")
    print(refseq.head())
    print()


    header("Filtering Bacteria Only...")
    refseq = refseq[(refseq["group"] == "bacteria") & (refseq["version_status"] == "latest") & (refseq["assembly_level"] == "Complete Genome")]
    refseq = refseq.loc[refseq.groupby(['assembly_accession', 'infraspecific_name'])['seq_rel_date'].idxmax()]
    #refseq = refseq.head()
    print(refseq.head())
    print()

    header("Adding Taxonomy Information...")
    taxa_db = gen_taxa_db(refseq, update)
    print(taxa_db.head())
    print()

    header(f"Saving to {outPath}...")

    taxa_db.to_csv(outPath, sep = "\t", index = False, mode = "w", float_format='%.0f')

    header("Done!")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("refseq_path")
    parser.add_argument("out_path")

    args = parser.parse_args()
    refseq_preprocessor(args.refseq_path, args.out_path)
