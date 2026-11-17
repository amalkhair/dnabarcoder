from Bio import SeqIO, Phylo
from pathlib import Path
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--raw", required=True)
parser.add_argument("--aligned", required=True)
parser.add_argument("--tree", required=True)
parser.add_argument("--marker", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

raw = {
    r.id: str(r.seq).upper()
    for r in SeqIO.parse(args.raw, "fasta")
}

aligned = {
    r.id: str(r.seq).upper()
    for r in SeqIO.parse(args.aligned, "fasta")
}

tree = Phylo.read(args.tree, "newick")

branches = {
    t.name: (t.branch_length or 0)
    for t in tree.get_terminals()
}

rows = []

for seq_id, seq in raw.items():

    aln = aligned[seq_id]

    rows.append({
        "tip": seq_id,
        "species": seq_id.split("|")[0].replace("_", " "),
        "raw_length": len(seq),
        "gap_percent": aln.count("-") / len(aln) * 100,
        "N_percent": seq.count("N") / len(seq) * 100,
        "branch_length": branches.get(seq_id, 0)
    })

df = pd.DataFrame(rows)

# Branch outliers
q1 = df["branch_length"].quantile(.25)
q3 = df["branch_length"].quantile(.75)
iqr = q3 - q1

branch_threshold = q3 + 3 * iqr

# Length outliers
lq1 = df["raw_length"].quantile(.25)
lq3 = df["raw_length"].quantile(.75)
liqr = lq3 - lq1

low_length = lq1 - 3 * liqr
high_length = lq3 + 3 * liqr

df["branch_outlier"] = (
    df["branch_length"] > branch_threshold
)

df["length_outlier"] = (
    (df["raw_length"] < low_length) |
    (df["raw_length"] > high_length)
)

df["high_N"] = df["N_percent"] > 5

df["needs_review"] = (
    df["branch_outlier"] |
    df["length_outlier"] |
    df["high_N"]
)

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(f"\n{args.marker} SEQUENCE QC")
print("=" * 80)

print(f"Sequences:             {len(df)}")
print(f"Median raw length:     {df['raw_length'].median():.0f}")
print(f"Length Q1-Q3:          {lq1:.0f}-{lq3:.0f}")
print(f"Length range allowed:  {low_length:.0f}-{high_length:.0f}")
print(f"Branch outliers:       {df['branch_outlier'].sum()}")
print(f"Length outliers:       {df['length_outlier'].sum()}")
print(f">5% N sequences:       {df['high_N'].sum()}")
print(f"Needs review:          {df['needs_review'].sum()}")

print("\nTOP LONG-BRANCH SEQUENCES\n")

print(
    df.sort_values("branch_length", ascending=False)[
        [
            "species",
            "raw_length",
            "gap_percent",
            "N_percent",
            "branch_length",
            "length_outlier"
        ]
    ]
    .head(20)
    .round(3)
    .to_string(index=False)
)

print(f"\nSaved: {out}")