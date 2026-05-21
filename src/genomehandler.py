from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from typing import List, Tuple
from iohandler import header
from tqdm import tqdm

def create_raw_multifasta(files: List[Path], out: Path):
    out.parent.mkdir(parents = True, exist_ok = True)
    with open(out, 'w', encoding='utf-8') as outfile:
        for file in tqdm(files, total=len(files)):
            with open(file, 'r', encoding='utf-8') as infile:
                for line in infile:
                    # Sanitization - Keeps only the protein accession in the header but also adds the file the protein came from
                    if line.startswith('>'):
                        protein_id = line[1:].split()[0]
                        genome_id = ".".join(file.name.split(".")[:-1])
                        outfile.write(f">{protein_id}|{genome_id}\n")
                    else:
                        outfile.write(line)
                
def create_diamond_db(infile: Path, outfile: Path):
    outfile.parent.mkdir(parents = True, exist_ok = True)
    header(f"Creating {outfile}.dmnd from {infile}")
    subprocess.run([
        "diamond",
        "makedb",
        "--in",
        str(infile),
        "--db",
        str(outfile)
    ])

def run_diamond(
    query: Path,
    db: Path,
    out: Path,
    threads: int = 4,
    max_target_seqs: int | None = None,
    coverage_thresh: float | None = None,
    identity_thresh: float | None = None,
    evalue: float | None = None,
):
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "diamond", "blastp",
        "-q", str(query),
        "-d", str(db),
        "-o", str(out),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "qlen", "qstart", "qend", "evalue", "bitscore",
        "--threads", str(threads),
    ]
    if max_target_seqs is not None:
        cmd.extend(["--max-target-seqs", str(max_target_seqs)])
    if coverage_thresh is not None:
        cmd.extend(["--query-cover", str(coverage_thresh)])
    if identity_thresh is not None:
        cmd.extend(["--id", str(identity_thresh)])
    if evalue is not None:
        cmd.extend(["--evalue", str(evalue)])

    subprocess.run(cmd, check=True)

def compute_pdt(
    tsv: Path,
    out: Path,
    query_ids: List[str],
    expected_genomes: List[str],
    identity_thresh: float,
    coverage_thresh: float,
):
    expected_genome_set = set(expected_genomes)
    total_genomes = len(expected_genome_set)

    # Start from full query universe so zero-hit genes are preserved
    pdt_df = pd.DataFrame({"qseqid": query_ids}).drop_duplicates()
    pdt_df["genomes_with_gene"] = 0
    pdt_df["total_genomes"] = total_genomes
    pdt_df["PDT"] = 0.0

    if not tsv.exists() or tsv.stat().st_size == 0:
        pdt_df.to_csv(out, sep="\t", index=False)
        return

    try:
        hits = pd.read_csv(
            tsv,
            sep="\t",
            header=None,
            names=["qseqid", "sseqid", "pident", "qlen", "qstart", "qend", "evalue", "bitscore"],
        )
    except pd.errors.EmptyDataError:
        pdt_df.to_csv(out, sep="\t", index=False)
        return

    if hits.empty:
        pdt_df.to_csv(out, sep="\t", index=False)
        return

    # Enforce Perl-like filtering explicitly
    hits["qcov"] = (hits["qend"] - hits["qstart"]).abs() + 1
    hits["qcov"] = 100.0 * hits["qcov"] / hits["qlen"].replace(0, pd.NA)
    hits = hits[(hits["pident"] >= identity_thresh) & (hits["qcov"] >= coverage_thresh)]

    if hits.empty:
        pdt_df.to_csv(out, sep="\t", index=False)
        return

    hits["genome"] = hits["sseqid"].astype(str).str.rsplit("|", n=1).str[-1]
    hits = hits[hits["genome"].isin(expected_genome_set)]

    if hits.empty:
        pdt_df.to_csv(out, sep="\t", index=False)
        return

    present = (
        hits[["qseqid", "genome"]]
        .drop_duplicates()
        .groupby("qseqid", as_index=False)
        .size()
        .rename(columns={"size": "genomes_with_gene"})
    )

    pdt_df = pdt_df.drop(columns=["genomes_with_gene"]).merge(present, on="qseqid", how="left")
    pdt_df["genomes_with_gene"] = pdt_df["genomes_with_gene"].fillna(0).astype(int)

    if total_genomes > 0:
        pdt_df["PDT"] = 100.0 * pdt_df["genomes_with_gene"] / total_genomes
    else:
        pdt_df["PDT"] = 0.0

    pdt_df.to_csv(out, sep="\t", index=False)


def find_hgt(row, pdt_thresh) -> Tuple[str, str]:
    """
    Returns (type, mode):
      type in {HGT, NATIVE, AMBIGUOUS}
      mode in {Recent, Ancient, Native}
    Mirrors Perl decision tree.
    """
    s_low, s_high = pdt_thresh["species"]["low"], pdt_thresh["species"]["high"]
    g_low, g_high = pdt_thresh["genus"]["low"], pdt_thresh["genus"]["high"]
    f_low, f_high = pdt_thresh["family"]["low"], pdt_thresh["family"]["high"]

    # species fallback column if present (Perl FLAG_pdt behavior)
    s_val = getattr(row, "PDT_species_for_rule", row.PDT_species_t2)

    g_t1, f_t1 = row.PDT_genus_t1, row.PDT_family_t1
    g_t2, f_t2 = row.PDT_genus_t2, row.PDT_family_t2

    if s_val <= s_low:
        return "HGT", "Recent"

    elif s_val >= s_high:
        if g_t2 >= g_high:
            if f_t2 >= f_high:
                return "NATIVE", "Native"
            elif f_t2 <= f_low:
                return "HGT", "Ancient"
            else:
                return "AMBIGUOUS", "Native"
        elif g_t2 <= g_high:
            if g_t2 >= g_low:
                if f_t1 <= f_low:
                    return "HGT", "Ancient"
                elif f_t1 >= f_high:
                    return "NATIVE", "Native"
                else:
                    return "AMBIGUOUS", "Native"
            elif g_t2 < g_low:
                return "HGT", "Ancient"
            else:
                return "AMBIGUOUS", "Native"
        else:
            return "AMBIGUOUS", "Native"

    elif (s_val < s_high) and (s_val > s_low):
        if g_t1 <= g_low:
            return "HGT", "Ancient"
        elif g_t1 >= g_high:
            return "NATIVE", "Native"
        elif (g_t1 > g_low) and (g_t1 < g_high):
            if f_t1 <= f_low:
                return "HGT", "Ancient"
            elif f_t1 >= f_high:
                return "NATIVE", "Native"
            else:
                return "AMBIGUOUS", "Native"
        else:
            return "AMBIGUOUS", "Native"

    return "AMBIGUOUS", "Native"