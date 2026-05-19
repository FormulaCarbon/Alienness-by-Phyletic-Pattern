from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
from itertools import batched
from iohandler import header, error
                
def download_dataset(accessions_file: Path, dehydrated_file: Path, db_path: Path):
    header("Downloading Dehydrated Dataset...")
    subprocess.run([
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        accessions_file,
        "--dehydrated",
        "--filename",
        dehydrated_file
        ]
    )
    print()
    
    header("Unzipping Dehydrated Dataset...")
    subprocess.run([
        "tar",
        "-xf",
        dehydrated_file,
        "-C",
        db_path,
        ], 
        shell = True
    )
    print()
    
    header("Rehydrating Dataset...")
    subprocess.run([
        "datasets",
        "rehydrate",
        "--directory",
        db_path,
        ]
    )
    print()