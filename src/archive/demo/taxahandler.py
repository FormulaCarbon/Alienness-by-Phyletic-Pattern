import pandas as pd
from ete3 import NCBITaxa

def add_families(target: pd.DataFrame) -> pd.DataFrame:
    ncbi = NCBITaxa()
    
    taxIDs = target["species_taxid"].dropna().astype(int).tolist()
    
    lineages = [ncbi.get_lineage(taxID) for taxID in taxIDs]
    families = []
    for lineage in lineages:
        familyIDs = []
        ranks = ncbi.get_rank(lineage)
        familyIDs = [familyID for familyID, rank in ranks.items() if rank == "family"]
        if familyIDs:
            families.append(ncbi.get_taxid_translator(familyIDs)[familyIDs[0]])
        else:
            families.append(None)
    target['family'] = families
    
    return target

def get_family(sciName: str) -> str | None:
    ncbi = NCBITaxa()
    
    taxID = ncbi.get_name_translator([sciName])[sciName][0]
    
    lineage = ncbi.get_lineage(taxID)
    ranks = ncbi.get_rank(lineage)
    familyIDs = [familyID for familyID, rank in ranks.items() if rank == "family"]
    if familyIDs:
        return ncbi.get_taxid_translator(familyIDs)[familyIDs[0]]
    return None


    

def search_organisms(refseq: pd.DataFrame, query: str, n: int = -1, exclude: str = "") -> pd.DataFrame:
    mask = refseq["organism_name"].str.contains(query, case = False)
    
    if exclude != "":
        mask &= ~refseq["organism_name"].str.contains(exclude, case = False)
    
    
    outDF = refseq.loc[mask]
    if n != -1:
        outDF = outDF.head(n)
    
    return outDF
    
def search_organisms_by_family(refseq: pd.DataFrame, query: str, n: int = -1, exclude: str = "") -> pd.DataFrame:
    mask = refseq["family"].str.contains(query, case = False)
    
    if exclude != "":
        mask &= ~refseq["organism_name"].str.contains(exclude, case = False)
    
    
    outDF = refseq.loc[mask]
    if n != -1:
        outDF = outDF.head(n)
    
    return outDF
    