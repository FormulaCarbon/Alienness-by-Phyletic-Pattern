from kmer_2level import read_fasta, gen_kmers, gen_proteome_kmers, presence_score
import argparse
from pathlib import Path
import pandas as pd

def classify_recent(speciesFrac: float, genusFrac: float, familyFrac: float, speciesThresh: float = 0.05, genusThresh: float = 0.3, familyThresh: float = 0.6):
    return (speciesFrac <= speciesThresh) and (genusFrac <= genusThresh) and (familyFrac >= familyThresh)

def classify_ancient(speciesFrac: float, genusFrac: float, familyFrac: float, speciesThresh: float = 0.6, genusThresh: float = 0.6, familyThresh: float = 0.6):
    return (speciesFrac >= speciesThresh) and (genusFrac >= genusThresh) and (familyFrac >= familyThresh)
        
def main(target: Path, bucketsPath: Path, family: str, genus: str, species: str, k: int, presenceThresh: float = 0.2) -> dict:
    speciesDir = bucketsPath / species
    genusDir = bucketsPath / genus
    familyDir = bucketsPath / family
    
    speciesFastas = sorted(speciesDir.glob("*.faa"))
    genusFastas = sorted(genusDir.glob("*.faa"))
    familyFastas = sorted(familyDir.glob("*.faa"))
    
    print("Generating Kmers")
    speciesIdx = {p.stem: gen_proteome_kmers(p, k) for p in speciesFastas}
    genusIdx = {p.stem: gen_proteome_kmers(p, k) for p in genusFastas}
    familyIdx = {p.stem: gen_proteome_kmers(p, k) for p in familyFastas}
    
    records = read_fasta(target)
    
    print("calculating scores")
    
    out = {}    
    
    for gene, seq in records.items():
        seqKmers = gen_kmers(seq, k)
        scores = {}
        present = [0, 0, 0]
        for name, genome_kmers in speciesIdx.items():
            score = presence_score(seqKmers, genome_kmers)
            scores[name] = score
            if score >= presenceThresh:
                present[0] += 1
                
        for name, genome_kmers in genusIdx.items():
            score = presence_score(seqKmers, genome_kmers)
            scores[name] = score
            if score >= presenceThresh:
                present[1] += 1
                
        for name, genome_kmers in familyIdx.items():
            score = presence_score(seqKmers, genome_kmers)
            scores[name] = score
            if score >= presenceThresh:
                present[2] += 1
                
        scores["present"] = sum(present)
        scores["species_frac"] = present[0]/len(speciesIdx)
        scores["genus_frac"] = present[1]/len(genusIdx)
        scores["family_frac"] = present[2]/len(familyIdx)
        scores["recent"] = classify_recent(scores["species_frac"], scores["genus_frac"], scores["family_frac"])
        scores["ancient"] = classify_ancient(scores["species_frac"], scores["genus_frac"], scores["family_frac"])
        out[gene] = scores
    
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("targetPath")
    parser.add_argument("family")
    parser.add_argument("genus")
    parser.add_argument("species")
    parser.add_argument("bucketsPath")
    parser.add_argument("outPath")
    args = parser.parse_args()
    
    targetPath = Path(args.targetPath)
    bucketsPath = Path(args.bucketsPath)
    
    res = main(targetPath, bucketsPath, args.family, args.genus, args.species, 6)
    
    df = pd.DataFrame.from_dict(res, orient="index")
    df.index.name = "gene"
    df = df.reset_index()
    
    df.to_csv(Path(args.outPath), sep = "\t")