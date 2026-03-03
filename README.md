recreation of APP but in python

to setup environment:
1. `py -m venv .venv`
2. `pip install -r requirements.txt`

to get the reference files needed:
1. download https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt and place into `src/ref/`
2. run `py src/refseq_preprocessor.py src/ref/assembly_summary_refseq.txt src/ref/BACTERIAL_LINEAGES.txt`

to download fasta files into buckets:
- if just needing one set of buckets: `py src/dlhandler.py taxid src/ref/BACTERIAL_LINEAGES.txt src/buckets/ --strain strain`, add `-e` if you want exclusive (T2) buckets
- if wanting both types of buckets `py src/dlhander.py taxid src/ref/BACTERIAL_LINEAGES.txt src/buckes/ --strain strain -c`

to generate a newick graph:
- `py src/taxahandler.py src/ref/BACTERIAL_LINEAGES.txt src/graph/tree.nwk`
