from collections import Counter, defaultdict
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import List
import pandas as pd
from tqdm import tqdm

from dlhandler import download_url, unzip_fastas

def process_fasta(faa: Path | str) -> Path:
    out = f"{faa}.modified"
    with open(faa, 'r') as infile, open(out, 'w') as outfile:
        counter = 0
        for line in infile:
            if line.startswith('>'):
                if counter > 0:
                    outfile.write("\n")
                counter += 1
                line = line.replace("-->", '').replace("->",'')
                arr1 = line.split(">")
                arr2 = arr1[1].split("] [")
                flag = 0
                real_header = line
                for key in arr2:
                    key = key.strip()
                    if "protein_id" in key:
                        acc_no = key.replace("protein_id", '')
                        acc_no = acc_no.split("=")[1]
                        if acc_no:
                            new_header = f"{acc_no}:g{counter}"
                            outfile.write(f">{new_header}\n")
                        else:
                            raise KeyError(line)
                        flag = 1
                if flag == 0:
                    real_header = real_header.strip()
                    arr3 = real_header.split(">")
                    fname = arr3[1]
                    new_header = f"{fname}:g{counter}"
                    outfile.write(f">{new_header}\n")
            else:
                outfile.write(line.strip())
    return Path(out)

def all_level_blast(acc: str, query: dict, lineages: pd.DataFrame, genomes_dir: Path, working_dir: Path, sid):
    kingdom, phylum, clas, order, family, genus, species = lineages.loc[lineages["accession"] == acc, ["kingdom_id", "phylum_id", "class_id", "order_id", "family_id", "genus_id", "species_id"]].iloc[0] 
    if sid != 0:
        species = sid
    order_under_consideration = order
    family_under_consideration = family
    genus_under_consideration = genus
    species_under_consideration = species
    flag_pdt = False
    # --- Type 1 ---
    
    print(kingdom, phylum, clas, order, family, genus, species)
    
    # Species

    s_matches = lineages.loc[(
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family))
        & (lineages["genus_id"].astype(str) == str(genus))
        & (lineages["species_id"].astype(str) == str(species_under_consideration))
    ), ["accession", "ftp_path"]]
    print(lineages)
    species_accs = s_matches["accession"].dropna().drop_duplicates().tolist()
    print(len(species_accs))
    #species_files = [genomes_dir / f"{acc}.faa" for acc in species_accs]
    species_files = []
    
    for acc in species_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            species_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 's1' / f"{acc}.faa.gz") == 0:
                species_files.append(working_dir / "temp" / 's1' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 's1' )
    
    if len(species_files) < 1:
        print("error, s t1 no files")
    else:
        print(f"species t1: {len(species_files)}")
        download_links_blastp_all(
            query, 
            species_files,
            len(species_files),
            "species",
            1,
            working_dir
        )
    
    # Genus
    g_matches = lineages.loc[(
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family))
        & (lineages["genus_id"].astype(str) == str(genus_under_consideration))
    ), ["accession", "ftp_path"]]
    
    genus_accs = g_matches["accession"].dropna().drop_duplicates().tolist()
    genus_files = []
    
    for acc in genus_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            genus_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 'g1' / f"{acc}.faa.gz") == 0:
                genus_files.append(working_dir / "temp" / 'g1' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 'g1' )
    
    if len(genus_files) < 1:
        print("error, g t1 no files")
    else:
        print(f"genus t1: {len(genus_files)}")
        download_links_blastp_all(
            query, 
            genus_files,
            len(genus_files),
            "genus",
            1,
            working_dir
        )
    
    # Family
    f_matches = lineages.loc[
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family_under_consideration)),
        ["accession", "ftp_path"],
    ]

    family_accs = f_matches["accession"].dropna().drop_duplicates().tolist()
    family_files = []
    
    for acc in family_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            family_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 'f1' / f"{acc}.faa.gz") == 0:
                family_files.append(working_dir / "temp" / 'f1' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 'f1' )

    if len(family_files) < 1:
        print("No family database found; hence going to order level (Type 1 inclusion)")

        # Order fallback
        o_matches = lineages.loc[
            (lineages["kingdom_id"].astype(str) == str(kingdom))
            & (lineages["phylum_id"].astype(str) == str(phylum))
            & (lineages["class_id"].astype(str) == str(clas))
            & (lineages["order_id"].astype(str) == str(order_under_consideration)),
            ["accession", "ftp_path"],
        ]
        order_accs = o_matches["accession"].dropna().drop_duplicates().tolist()
        order_files = [genomes_dir / f"{acc}.faa" for acc in order_accs if (genomes_dir / f"{acc}.faa").exists()]

        if len(order_files) < 1:
            print("error, order t1 no files")
        else:

            download_links_blastp_all(
                query,
                order_files,
                len(order_files),
                "order",
                1,
                working_dir,
            )
    else:
        print(f"family t1: {len(family_files)}")
        download_links_blastp_all(
            query,
            family_files,
            len(family_files),
            "family",
            1,
            working_dir,
        )
    
    # --- Type 2 ---
    
    # Species
    s_matches = lineages.loc[
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family))
        & (lineages["genus_id"].astype(str) == str(genus))
        & (lineages["species_id"].astype(str) == str(species_under_consideration))
        & (lineages["strain"].astype(str) != str(query["strain"])),
        ["accession", "ftp_path"]]
    
    species_accs = s_matches["accession"].dropna().drop_duplicates().tolist()
    species_files = []
    
    for acc in species_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            species_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 's2' / f"{acc}.faa.gz") == 0:
                species_files.append(working_dir / "temp" / 's2' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 's2' )
    
    if len(species_files) < 1:
        print("error, s t2 no files")
        flag_pdt = True
    else:
        print(f"species t2: {len(species_files)}")
        download_links_blastp_all(
            query, 
            species_files,
            len(species_files),
            "species",
            2,
            working_dir
        )

    print(kingdom, phylum, clas, order, family, genus_under_consideration, species_under_consideration, query['tax_id'])
    # Genus
    g_matches = lineages.loc[
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family))
        & (lineages["genus_id"].astype(str) == str(genus_under_consideration))
        & (lineages["species_id"].astype(str) != str(species_under_consideration)),
        #& (lineages["tax_id"].astype(str) != str(query["tax_id"])),
        ["accession", "ftp_path"]]
    print(len(g_matches))

    
    genus_accs = g_matches["accession"].dropna().drop_duplicates().tolist()
    print(len(genus_accs))
    genus_files = []
    
    for acc in genus_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            genus_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 'g2' / f"{acc}.faa.gz") == 0:
                genus_files.append(working_dir / "temp" / 'g2' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 'g2' )
    
    if len(genus_files) < 1:
        print("error, g t2 no files")
    else:
        print(f"genus t2: {len(genus_files)}")
        download_links_blastp_all(
            query, 
            genus_files,
            len(genus_files),
            "genus",
            2,
            working_dir
        )
    
    # Family
    f_matches = lineages.loc[
        (lineages["kingdom_id"].astype(str) == str(kingdom))
        & (lineages["phylum_id"].astype(str) == str(phylum))
        & (lineages["class_id"].astype(str) == str(clas))
        & (lineages["order_id"].astype(str) == str(order))
        & (lineages["family_id"].astype(str) == str(family_under_consideration))
        & (lineages["genus_id"].astype(str) != str(genus_under_consideration))
        & (lineages["species_id"].astype(str) != str(species_under_consideration)),
        #& (lineages["tax_id"].astype(str) != str(query["tax_id"])),
        ["accession", "ftp_path"],
    ]

    family_accs = f_matches["accession"].dropna().drop_duplicates().tolist()
    family_files = []
    for acc in family_accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            family_files.append(file)
        else:
            if download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / 'f2' / f"{acc}.faa.gz") == 0:
                family_files.append(working_dir / "temp" / 'f2' / f"{acc}.faa")
            else:
                print('error downloading ', acc)
                continue
    unzip_fastas(working_dir / "temp" / 'f2' )

    if len(family_files) < 1:
        print("No family database found; hence going to order level (Type 2 inclusion)")

        # Order fallback
        o_matches = lineages.loc[
            (lineages["kingdom_id"].astype(str) == str(kingdom))
            & (lineages["phylum_id"].astype(str) == str(phylum))
            & (lineages["class_id"].astype(str) == str(clas))
            & (lineages["order_id"].astype(str) == str(order_under_consideration))
            & (lineages["family_id"].astype(str) != str(family_under_consideration))
            & (lineages["genus_id"].astype(str) != str(genus_under_consideration))
            & (lineages["species_id"].astype(str) != str(species_under_consideration)),
            #& (lineages["tax_id"].astype(str) != str(query["tax_id"])),
            ["accession", "ftp_path"],
        ]
        order_accs = o_matches["accession"].dropna().drop_duplicates().tolist()
        order_files = [genomes_dir / f"{acc}.faa" for acc in order_accs if (genomes_dir / f"{acc}.faa").exists()]

        if len(order_files) < 1:
            print("error, order t2 no files")
        else:

            download_links_blastp_all(
                query,
                order_files,
                len(order_files),
                "order",
                2,
                working_dir,
            )
    else:
        print(f"family t2: {len(family_files)}")
        download_links_blastp_all(
            query,
            family_files,
            len(family_files),
            "family",
            2,
            working_dir,
        )
    return flag_pdt
    
def download_links_blastp_all(query: dict, files: List[Path], size, level, t, working_dir: Path):
    fasta = working_dir / f"{level}_{t}_raw_database.fasta"
    map_file = working_dir / f"{level}_{t}_Accession_no_species.txt"
    print(f"creating multifasta {fasta}")
    create_raw_multifasta(files, fasta, map_file)
    
    print("creating blast db")
    subprocess.run([
        "makeblastdb",
        "-in", str(fasta),
        "-parse_seqids",
        "-dbtype", "prot"
    ], check=True)
    
    threads = detect_cpu_cores()
    print(f"running blast with {threads} threads")
    subprocess.run([
        "blastp",
        "-num_threads", str(threads),
        "-query", str(query["file_path"]),
        "-db", str(fasta),
        "-out", str(working_dir / f"{level}_{t}_output_blasthits"),
        "-outfmt", '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore',
        "-max_target_seqs", str(size)
    ], check=True)
    
def download_links_blastp_all_diamond(query: dict, files: List[Path], size, level, t, working_dir: Path):
    fasta = working_dir / f"{level}_{t}_raw_database.fasta"
    map_file = working_dir / f"{level}_{t}_Accession_no_species.txt"
    outfile = working_dir / f"{level}_{t}_db"
    print(f"creating db {fasta}")
    #create_raw_multifasta(files, fasta, map_file)
    
    print("creating diamond db")
    outfile.parent.mkdir(parents = True, exist_ok = True)
    print(f"Creating {outfile}.dmnd from {fasta}")
    subprocess.run([
        "diamond",
        "makedb",
        "--in",
        str(fasta),
        "--db",
        str(outfile),
        "--verbose"
    ], check= True)
    
    threads = detect_cpu_cores()
    print(f"running diamond with {threads} threads")
    cmd = [
        "diamond", "blastp",
        "-q",  str(query["file_path"]),
        "-d", str(outfile),
        "-o", str(working_dir / f"{level}_{t}_output_blasthits"),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "qcovhsp", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "--threads", str(threads),
        "--max-target-seqs", str(size),
        "--ultra-sensitive",
        "--comp-based-stats", "4"
    ]
    
    print(cmd)
    subprocess.run(cmd, check=True)

    
def download_links_blastp_all_mmseqs2(query: dict, files: List[Path], size, level, t, working_dir: Path):
    fasta = working_dir / f"{level}_{t}_raw_database.fasta"
    map_file = working_dir / f"{level}_{t}_Accession_no_species.txt"
    print(f"creating multifasta {fasta}")
    create_raw_multifasta(files, fasta, map_file)

    db_path = working_dir / f"{level}_{t}_db"
    query_db_path = working_dir / f"{level}_{t}_query_db"
    result_db_path = working_dir / f"{level}_{t}_result_db"
    output_path = working_dir / f"{level}_{t}_output_blasthits"
    tmp_path = working_dir / f"{level}_{t}_tmp"
    tmp_path.mkdir(exist_ok=True)

    print("creating mmseqs2 db")
    subprocess.run(["mmseqs", "createdb", str(fasta), str(db_path)], check=True)
    subprocess.run(["mmseqs", "createdb", str(query["file_path"]), str(query_db_path)], check=True)

    threads = detect_cpu_cores()
    print(f"running mmseqs2 with {threads} threads")
    subprocess.run([
        "mmseqs", "search",
        str(query_db_path), str(db_path), str(result_db_path), str(tmp_path),
        "--threads", str(threads),
        "--max-seqs", str(size),
        "-s", "7.5",
    ], check=True)

    subprocess.run([
        "mmseqs", "convertalis",
        str(query_db_path), str(db_path), str(result_db_path), str(output_path),
        "--threads", str(threads),
        "--format-mode", "0",
        "--format-output", "query,target,pident,qcov,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    ], check=True)
    
    shutil.rmtree(tmp_path, ignore_errors=True)
    
    

def create_raw_multifasta(files: List[Path], out: Path, map_file: Path):
    out.parent.mkdir(parents = True, exist_ok = True)
    with open(out, 'w', encoding='utf-8') as outfile, open(map_file, "w", encoding="utf-8") as m:
        counter = 1
        for file in tqdm(files, total=len(files)):
            filename = file.stem.replace('.', '_')
            with open(file, 'r', encoding='utf-8') as infile:
                for line in infile:
                    if line.startswith('>'):
                        protein_id = line[1:].split()[0]
                        outfile.write(f">{protein_id}:g{counter}\n")
                        m.write(f"{protein_id}:g{counter}\t{filename}\n")
                    else:
                        outfile.write(line)
            counter += 1

def detect_cpu_cores() -> int:
    # Mirror Perl logic as close as possible
    # If /proc/cpuinfo: grep -c ^processor
    if os.path.exists("/proc/cpuinfo"):
        try:
            out = subprocess.check_output("grep -c ^processor /proc/cpuinfo", shell=True, text=True)
            n = int(out.strip())
            if n > 0:
                return n
        except Exception:
            pass

    if sys.platform == "darwin":
        try:
            out = subprocess.check_output("sysctl -n hw.ncpu", shell=True, text=True)
            n = int(out.strip())
            if n > 0:
                return n
        except Exception:
            pass

    if os.name == "nt":
        v = os.environ.get("NUMBER_OF_PROCESSORS", "").strip()
        if v.isdigit():
            return int(v)

    print("Cannot determine the number of CPUs. Do single threading.\n")
    return 1

from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd

def output_blasthit_processor(
    t,
    hits: Path,
    map_file: Path,
    thresh_pident: int,
    level: str,
    query_file: Path,
    working_dir: Path,
    config
):
    names = "qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore".split()
    blasthits = pd.read_csv(hits, sep="\t", header=None, names=names)

    mappings = (
        pd.read_csv(map_file, sep="\t", header=None, names=["prot_acc", "gen_acc"])
          .dropna()
          .set_index("prot_acc")["gen_acc"]
          .to_dict()
    )

    with open(query_file, "r", encoding="utf-8", errors="replace") as f:
        sequences = f.readlines()

    working_dir.mkdir(parents=True, exist_ok=True)

    out1_path = working_dir / f"Main_blast_output_counts_pident_{thresh_pident}_{level}_{t}.txt"
    out2_path = working_dir / f"{thresh_pident}_calculations_{level}_{t}.txt"
    out3_path = working_dir / f"FOR_INPUT_{level}_{t}.txt"

    # Perl: uniq genomes from mapping file
    uniq_all_genomes = sorted(set(mappings.values()))

    # Perl: query gene IDs from fasta headers (everything after '>')
    JSCB_accno_all = []
    for s in sequences:
        s = s.strip()
        if s.startswith(">"):
            # safer: take first token (optional); if you want Perl-exact, keep s[1:]
            gene = s[1:]
            JSCB_accno_all.append(gene)

    # Perl: initialize all gene x genome combos to 0
    blasthits_counter = defaultdict(Counter)
    for gene in JSCB_accno_all:
        for genome in uniq_all_genomes:
            blasthits_counter[gene][genome] = 0

    # Perl: fill out hits as 1 if thresholds pass
    for hit in blasthits.itertuples(index=False):
        query = str(hit.qseqid)
        subject = str(hit.sseqid)
        pident = float(hit.pident)
        qcov = float(hit.qcovs)   # IMPORTANT: column is qcovs

        target = mappings.get(subject)
        if target is None:
            # subject id not found in mapping file
            continue

        if pident >= thresh_pident and qcov >= config['coverage_threshold']:
            blasthits_counter[query][target] = 1

    # --- write tables (OUTPUT1 and OUTPUT3 are identical in Perl) ---
    header_cols = ["genes"] + uniq_all_genomes

    with open(out1_path, "w", encoding="utf-8") as out1, \
         open(out3_path, "w", encoding="utf-8") as out3:

        

        out1.write("\t".join(header_cols) + "\n")
        out3.write("\t".join(header_cols) + "\n")

        n_genomes = len(uniq_all_genomes)

        for gene in sorted(blasthits_counter.keys()):
            row_vals = [str(blasthits_counter[gene][g]) for g in uniq_all_genomes]
            line = "\t".join([gene] + row_vals) + "\n"

            out1.write(line)
            out3.write(line)

    
    with open(out3_path, "r", encoding="utf-8") as infile, \
         open(out2_path, "w", encoding="utf-8") as out2:

        out2.write("Gene\tno. of genomes\tcount in which present\tfraction of genomes with the gene\n")
        
        for line in infile:
            line = line.strip()
            
            if re.match(r"^genes", line):
                continue
            
            arr = line.split('\t')
            gene = arr[0].strip()
            
            gene_counter = 0
            
            for i in range(0, len(arr) - 1):
                try:
                    if int(arr[i + 1]) == 1:
                        gene_counter += 1
                except ValueError:
                    pass
                    
            size = len(arr)-1
            frac = (float(gene_counter) / float(size)) * 100
            
            out2.write(f"{gene}\t{size}\t{gene_counter}\t{frac}\n")
            
def alien_genes_finder(
    flag, 
    working_dir: Path,
    calc_s1: Path,
    calc_g1: Path,
    calc_f1: Path,
    calc_s2: Path,
    calc_g2: Path,
    calc_f2: Path,
    query_file : Path,
    config
    ):
    
    output11 = working_dir / "phyletic_pattern_analysis_native_alien_genes_full_record.txt"
    output12 = working_dir / "APP_Alien_genes.txt"
    
    pdts_s1 = {}
    if flag == 1:
        with open(calc_s1, 'r', encoding="utf-8") as species_type1:
            for line in species_type1:
                line = line.strip()
                if re.match(r"^Gene", line):
                    continue
                arr = line.split("\t")
                acc = arr[0].strip()
                pdt = arr[3].strip()
                pdts_s1[acc] = {"pdt": float(pdt)}
    
    pdts_s2 = {}
    if flag == 0:
        with open(calc_s2, 'r', encoding="utf-8") as species_type2:
            for line in species_type2:
                line = line.strip()
                if re.match(r"^Gene", line):
                    continue
                arr = line.split("\t")
                acc = arr[0].strip()
                pdt = arr[3].strip()
                pdts_s2[acc] = {"pdt": float(pdt)}
                
    accession_nos = {}
    counter = 1
    with open(query_file, 'r') as q:
        for line in q:
            if line.startswith(">"):
                g = line.split(">")[1].strip()
                accession_nos[counter] = g
                counter += 1
    
    pdts_g1 = {}
    with open(calc_g1, 'r', encoding="utf-8") as genus_type1:
        for line in genus_type1:
            line = line.strip()
            if re.match(r"^Gene", line):
                continue
            arr = line.split("\t")
            acc = arr[0].strip()
            pdt = arr[3].strip()
            pdts_g1[acc] = {"pdt": float(pdt)}
            
    pdts_f1 = {}
    with open(calc_f1, 'r', encoding="utf-8") as family_type1:
        for line in family_type1:
            line = line.strip()
            if re.match(r"^Gene", line):
                continue
            arr = line.split("\t")
            acc = arr[0].strip()
            pdt = arr[3].strip()
            pdts_f1[acc] = {"pdt": float(pdt)}
    
    # Type 2
    
    pdts_g2 = {}
    with open(calc_g2, 'r', encoding="utf-8") as genus_type2:
        for line in genus_type2:
            line = line.strip()
            if re.match(r"^Gene", line):
                continue
            arr = line.split("\t")
            acc = arr[0].strip()
            pdt = arr[3].strip()
            pdts_g2[acc] = {"pdt": float(pdt)}
            
    pdts_f2 = {}
    with open(calc_f2, 'r', encoding="utf-8") as family_type2:
        for line in family_type2:
            line = line.strip()
            if re.match(r"^Gene", line):
                continue
            arr = line.split("\t")
            acc = arr[0].strip()
            pdt = arr[3].strip()
            pdts_f2[acc] = {"pdt": float(pdt)}
            
    out = {}
    species_data = pdts_s1 if flag == 1 else pdts_s2
    
    for acc in species_data:
        pdt_species_t1 = float(pdts_s1.get(acc, {}).get("pdt", 0.0))
        pdt_genus_t1 = float(pdts_g1.get(acc, {}).get("pdt", 0.0))
        pdt_family_t1 = float(pdts_f1.get(acc, {}).get("pdt", 0.0))

        # exclusion
        pdt_species_t2 = float(pdts_s2.get(acc, {}).get("pdt", 0.0))
        pdt_genus_t2 = float(pdts_g2.get(acc, {}).get("pdt", 0.0))
        pdt_family_t2 = float(pdts_f2.get(acc, {}).get("pdt", 0.0))
        
        
        ##Decision rules; 
        #"type" => Native = 1; Alien = 0
        #"mode" => Ancient = 1; Recent = 0; Native = 2; Native values assigned to avoid error warning
        #"PDT_species_T1_or_T2" => When exlusion data is not available, type 1 data is used at species level
        # TODO: add config support
        if flag == 1:
            print("using incl data for s")
            pdt_species = pdt_species_t1
        elif flag == 0:
            print("using excl data for s")
            pdt_species = pdt_species_t2
            
        if pdt_species <= 30:
            out[acc] = {"type": 0, "mode": 0}
        
        elif pdt_species >= 80:
            if pdt_genus_t2 >= 70:
                if pdt_family_t2 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t2 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t2 > 30 and pdt_family_t2 < 40:
                    out[acc] = {"type": "Ambiguous1", "mode": 2}
                else:
                    out[acc] = {"type": "Error1"}
            elif pdt_genus_t2 <= 70 and pdt_family_t2 >= 30:
                if pdt_family_t1 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t1 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t1 > 30 and pdt_family_t1 < 40:
                    out[acc] = {"type": "Ambiguous2", "mode": 2}
                else:
                    out[acc] = {"type": "Error2"}
            elif pdt_genus_t2 < 30:
                out[acc] = {"type": 0, "mode": 1}
            else:
                out[acc] = {"type": "Error3"}
        
        elif pdt_species < 80 and pdt_species > 30:
            if pdt_genus_t1 <= 30:
                out[acc] = {"type": 0, "mode": 1}
            elif pdt_genus_t1 >= 70:
                out[acc] = {"type": 1, "mode": 2}
            elif pdt_genus_t1 > 30 and pdt_genus_t1 < 70:
                if pdt_family_t1 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t1 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t1 > 30 and pdt_family_t1 < 40:
                    out[acc] = {"type": "Ambiguous3", "mode": 2}
                else:
                    out[acc] = {"type": "Error4"}
        
        else:
            out[acc] = {"type": "Error"}
        with open(output11, 'w') as o11, open(output12, 'w') as o12:
            for a in sorted(out.keys()):
                b = pdts_s1.get(a, {}).get("pdt", "")
                c = pdts_g1.get(a, {}).get("pdt", "")
                d = pdts_f1.get(a, {}).get("pdt", "")
                e = pdts_s2.get(a, {}).get("pdt", "")
                f = pdts_g2.get(a, {}).get("pdt", "")
                g = pdts_f2.get(a, {}).get("pdt", "")
                h = out[a].get("type", "")
                t = out[a].get("mode", "")
                o11.write(f"{a}\t{b}\t{c}\t{d}\t{e}\t{f}\t{g}\t{h}\t{t}\n")
                if h == 0:
                    if t == 1:
                        o12.write(f"{a}\t{t}\tAncient\n")
                    elif t == 0:
                        o12.write(f"{a}\t{t}\tRecent\n")

def marker_gene_enrichment(query_file, output_file):

    with open(query_file, 'r') as q:
        all_genes = q.readlines()
        
    with open(output_file, 'r') as o:
        alien_genes = o.readlines()
        
    total_genes_on_genome = 0
    all_gene_enrichment = []
    for i in range(0, len(all_genes) - 1, 2):
        line = all_genes[i]
        all_gene_enrichment.append(line.split(">")[1].strip())
        total_genes_on_genome += 1
        
    total_alien_genes = 0
    alien_gene_enrichment = []
    
    for g in alien_genes:
        g = g.strip()
        arr = g.split("\t")
        gene = arr[0]
        alien_gene_enrichment.append(gene)
        total_alien_genes += 1
        
    print(total_alien_genes)
    print(total_genes_on_genome)
    
    # implement rest later