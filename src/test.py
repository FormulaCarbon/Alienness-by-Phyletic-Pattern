import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm

lineages = pd.read_csv(Path("H:/GENOME_DB/BACTERIAL_LINEAGES.txt"), sep = '\t')

fasta_dir = Path("H:/GENOME_DB/GENOMES") / "ncbi_dataset" / "data"
file_count = sum(1 for x in fasta_dir.rglob('*') if x.is_file()) - 2


with open('log.txt', 'w') as f:
        for accession in tqdm(fasta_dir.iterdir(), total=file_count):
            try:
                clean = accession.name
                row = list(lineages.loc[lineages["accession"] == clean, ["order", "family", "genus", "species"]].iloc[0])
                filename = os.listdir(accession)[0]
                outpath = Path("H:/GENOME_DB/GENOMES") / row[0] / row[1] / row[2] / row[3] / filename 
                source = accession / filename
                outpath.parent.mkdir(exist_ok = True, parents = True)
                source.rename(outpath)
                os.removedirs(accession)
            except Exception as e:
                f.write(str(e) + "\n")
        f.close()  
