from dlhandler import *

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("queryPath")
parser.add_argument("dbPath")
parser.add_argument("outPath")
parser.add_argument("--threads", type=int, default=4, help="Number of threads for DIAMOND")
parser.add_argument("--total_genomes", type=int, default=None, help="Total genomes in the DB for PDT calculation")
parser.add_argument("--identity", type=float, default=50.0, help="Minimum percent identity")

args = parser.parse_args()

queryPath = Path(args.queryPath)
dbPath = Path(args.dbPath)
outPath = Path(args.outPath)

print(f"Running DIAMOND: {queryPath} vs {dbPath}")
run_diamond(queryPath, dbPath, outPath, threads=args.threads)

print("Computing PDT...")
pdt_dict = compute_pdt(outPath, args.identity, 0.7, args.total_genomes)

# save PDT to CSV
pdt_df = pd.DataFrame(list(pdt_dict.items()), columns=["gene", "PDT"])
pdt_csv = outPath.with_suffix(".pdt.csv")
pdt_df.to_csv(pdt_csv, index=False)
print(f"PDT saved to {pdt_csv}")
