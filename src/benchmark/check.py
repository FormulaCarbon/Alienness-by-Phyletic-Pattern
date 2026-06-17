import argparse
from pathlib import Path
import pandas as pd
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("outputs")
args = parser.parse_args()

outputs = Path(args.outputs)
alien_pred = []

for item in outputs.iterdir():
    if item.is_dir() and item.name[:3] == "GCF":
        alien_genes_df = pd.read_csv(item / "APP_Alien_genes.txt", sep = '\t', header=None, names = ["gene", "type", "age"])
        genes = alien_genes_df['gene'].apply(lambda x: f"{x.split(':')[0]}|{item.name}").to_list()
        alien_pred += [gene for gene in genes if gene.startswith('WP')]

with open("src/benchmark/detected.txt", 'w') as f:
    f.write('\n'.join(alien_pred))
"""
alien_true = []
with open("src/benchmark/alien_genes.faa", 'r') as f:
    for line in f:
        if line.startswith('>WP'):
            alien_true.append(line[1:])
with open("src/benchmark/alien.txt", 'w') as f:
    f.writelines(alien_true)
"""
alien_true = []
with open("src/benchmark/alien.txt", 'r') as f:
    for line in f:
        alien_true.append(line.strip())
        
print(alien_true[:5], len(alien_true), len(set(alien_true)))
        
alien_pred = []
with open("src/benchmark/detected.txt", 'r') as f:
    for line in f:
        alien_pred.append(line.strip())
print(alien_pred[:5], len(alien_pred), len(set(alien_pred)))
not_found = list(set(alien_true) - set(alien_pred))
print("not found:", len(not_found))
            
incorrect = list(set(alien_pred) - set(alien_true))
print("incorrect:", len(incorrect))

correct = list(set(alien_true) & set(alien_pred))
print("correct:", len(correct))