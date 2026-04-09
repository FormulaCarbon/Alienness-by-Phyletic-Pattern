# Alieness by Phyletic Pattern - Python Edition

recreation of APP but in python

make sure python is installed. Has been tested on 3.14. depending on the install `py` may need to be replaced by `python`

## Setup
1. `py -m venv .venv`
2. `.venv\scripts\activate`
3. `pip install -r requirements.txt`

## Execution
### Creating Reference Files
For bacterial lineages:
1. download https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
2. run `py src/refseq_preprocessor.py path/to/downloaded/file src/ref/BACTERIAL_LINEAGES.txt`
This will create a tab-seperated values file containing the full taxonomic ranks for all latest complete bacterial sequences present in NCBI's refseq (~54,000), as well as there taxonomic ID (taxid), rank IDs, and the download link for their FASTA.

For IslandPick:
- Wherever BACTERIAL_LINEAGES.txt is used in the following snippets, replace it with the path to `islandpick_lineages.txt` in `src/examples/`
Instructions on creating this file will be provided soon, for now just use the given one. Note that this does not include all items from islandpick, as it gets the data of the bacteria by referencing taxids with the NCBI database, so any bacteria in IslandPick but not NCBI's refseq will not be in this file.

### Downloading relatives
1. `py src/dlhander.py target_taxid src/ref/BACTERIAL_LINEAGES.txt src/buckets/ -c`, use `--strain "strain name"` if the target bacteria is a strain, use `--speciesdist INT`, `--genusdist int`, and `--familydist int` to set a max on how many relative are downloaded. Set to 0 to download all.
This will download two buckets, T1 and T2 as described in the paper. This will also create the databases for diamond

### Running DIAMOND:
1. `py src/APP.py target_fasta.faa database.dmnd output.tsv --identity INT --coverage FLOAT --total_genomes INT`
the total genomes will be equivalent to the number of genomes in the corresponding folder of the database.

### Graphing
to generate a newick graph:
- `py src/taxahandler.py src/ref/BACTERIAL_LINEAGES.txt src/graphs/tree.nwk`
- 
## DEV STUFF
get taxids from accessions
`cat acc.txt | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > taxids.txt`
