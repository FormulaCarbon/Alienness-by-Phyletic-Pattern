import subprocess
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument('regions_data', type = str)
parser.add_argument('alien_genes')
parser.add_argument('genome_db')
parser.add_argument('out_dir')

args = parser.parse_args()
out_dir = Path(args.out_dir)

regions_data = {}
with open(args.regions_data, 'r') as f:
    regions_data = json.load(f)
    
accs = [regions_data[name]['accession'] for name in regions_data]

for acc in accs:
    cmd = [
        'py',
        'src/app2.py',
        str(acc),
        str(args.genome_db),
        '--dir', str(out_dir / acc)
    ]
    
    subprocess.run(cmd, check=True)
    
    