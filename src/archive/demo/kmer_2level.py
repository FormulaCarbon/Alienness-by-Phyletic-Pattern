from Bio import SeqIO
from pathlib import Path
import pandas as pd

import argparse

def read_fasta(filePath: Path) -> dict:
    records = {}
    return {record.id: str(record.seq) for record in SeqIO.parse(filePath, "fasta")}

def gen_kmers(seq: str, k: int) -> set[str]:
    return set() if len(seq) < k else {seq[i:i+k] for i in range(len(seq)-k+1)}

def gen_proteome_kmers(fastaPath: Path, k: int) -> set[str]:
    out = set()
    for seqID, seq in read_fasta(fastaPath).items():
        out.update(gen_kmers(seq, k))
    return out
    
def presence_score(seqKmers: set[str], kmers: set[str]):
    return 0.0 if not seqKmers else len(seqKmers & kmers) / len(seqKmers)

def main(target: Path, bucketsPath: Path, genus: str, species: str, k: int, presenceThresh: float = 0.2) -> dict:
    speciesDir = bucketsPath / species
    genusDir = bucketsPath / genus
    
    speciesFastas = sorted(speciesDir.glob("*.faa"))
    genusFastas = sorted(genusDir.glob("*.faa"))
    
    print("Generating Kmers")
    speciesIdx = {p.stem: gen_proteome_kmers(p, k) for p in speciesFastas}
    genusIdx = {p.stem: gen_proteome_kmers(p, k) for p in genusFastas}
    
    records = read_fasta(target)
    
    print("calculating scores")
    
    out = {}    
    
    for gene, seq in records.items():
        seqKmers = gen_kmers(seq, k)
        scores = {}
        present = 0
        for name, genome_kmers in speciesIdx.items():
            score = presence_score(seqKmers, genome_kmers)
            scores[name] = score
            if score >= presenceThresh:
                present += 1
                
        for name, genome_kmers in genusIdx.items():
            score = presence_score(seqKmers, genome_kmers)
            scores[name] = score
            if score >= presenceThresh:
                present += 1
                
        scores["present"] = present
        out[gene] = scores
    
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("targetPath")
    parser.add_argument("genus")
    parser.add_argument("species")
    parser.add_argument("bucketsPath")
    parser.add_argument("outPath")
    args = parser.parse_args()
    
    targetPath = Path(args.targetPath)
    bucketsPath = Path(args.bucketsPath)
    
    res = main(targetPath, bucketsPath, args.genus, args.species, 6)
    
    df = pd.DataFrame.from_dict(res, orient="index")
    df.index.name = "gene"
    df = df.reset_index()
    
    df.to_csv(Path(args.outPath), sep = "\t")
