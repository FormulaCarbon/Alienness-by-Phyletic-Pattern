import yaml
from pathlib import Path

def create_config(db_path: Path):
    data = {
        "last_updated": "0001-01-01",
        "genomes_path": "GENOMES/",
        "lineage_path": "BACTERIAL_LINEAGES.txt",
        "refseq_path": "assembly_summary_refseq.txt",
        "log_path": "log.txt"
    }
    with open(db_path / "config.yaml", 'w') as f:
        yaml.dump(data, f, default_flow_style = False)