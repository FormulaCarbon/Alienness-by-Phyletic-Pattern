from pathlib import Path
import subprocess
import pandas as pd

def make_diamond_db(genomeFile: Path, outPath: Path):
    subprocess.run([
        "diamond", "makedb",
        "--in", str(genomeFile),
        "-d", str(outPath)
    ], check=True)
    
def run_diamond(query: Path, db: Path, out: Path, threads: int = 4): # chatgpt, redo it
    subprocess.run([
        "diamond", "blastp",
        "-q", str(query),
        "-d", str(db),
        "-o", str(out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen",
        "--threads", str(threads),
        "--fast"
    ], check=True)
    
def compute_pdt(tsv: Path, identity: float, coverage: float, total_genomes: int): # chatgpt, redo it
    df = pd.read_csv(tsv, sep="\t", header=None,
        names=["qseqid","sseqid","pident","length","qlen"])

    # coverage
    df["qcov"] = df["length"] / df["qlen"]

    # filter like APP
    df = df[(df["pident"] >= identity) & (df["qcov"] >= coverage)]

    # extract genome
    df["genome"] = df["sseqid"].str.split("|").str[0]

    # compute PDT
    pdt = (
        df.groupby("qseqid")["genome"]
        .nunique()
        / total_genomes
    )

    return pdt.to_dict()