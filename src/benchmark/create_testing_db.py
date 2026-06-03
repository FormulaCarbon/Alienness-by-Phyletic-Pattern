import pandas as pd
import json

with open('src/benchmark/accessions.txt', 'r') as f:
    accs = json.load(f)

def format_coordinates(coordinates:str):
    l = coordinates.split()
    out = []
    for c in l:
        s = c.split('..')
        out.append(tuple(map(lambda x: int(x), s)))
    return out

def nc_to_gcf(ncs: list):
    out = []
    for nc in ncs:
        out.append(accs[nc])
    return out

data = pd.read_csv("src/benchmark/18_db.csv", header=1)
data = data.drop(data.columns[0], axis = 1).drop(columns="Reference").drop(columns="Evidence")
data = data.rename(columns = {'Genomes':'genomes', 'Accession no':'accession', 'Genomic Island co-ordinates':'coordinates'})
data['coordinates'] = data['coordinates'].apply(format_coordinates)
data['accession'] = data['accession'].apply(lambda x: x.split(", ")).apply(nc_to_gcf)
out = {}
for line in data.itertuples():
    out[line.genomes] = {"accession": line.accession[0], "coordinates": line.coordinates}
with open('src/benchmark/data.json', 'w') as outfile:
    json.dump(out, outfile)