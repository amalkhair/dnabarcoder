from Bio import Phylo
from pathlib import Path
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tree", required=True)
parser.add_argument("--marker", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

tree_path = Path(args.tree)
out_dir = Path(args.out)
out_dir.mkdir(parents=True, exist_ok=True)

tree = Phylo.read(tree_path, "newick")

# --------------------------------------------------
# Parse species from our FASTA headers
# Example:
# Fusarium_oxysporum|ITS|rep1|ABC123
# --------------------------------------------------

tips = tree.get_terminals()

def get_species(tip):
    name = tip.name.split("|")[0]
    return name.replace("_", " ")

species_to_tips = {}

for tip in tips:
    species = get_species(tip)
    species_to_tips.setdefault(species, []).append(tip.name)

all_tip_names = {tip.name for tip in tips}

# --------------------------------------------------
# Build every bipartition in the unrooted tree
# --------------------------------------------------

splits = set()

for clade in tree.find_clades(order="preorder"):

    if clade == tree.root:
        continue

    side = {
        tip.name
        for tip in clade.get_terminals()
    }

    other = all_tip_names - side

    splits.add(frozenset(side))
    splits.add(frozenset(other))

# --------------------------------------------------
# Species monophyly
# --------------------------------------------------

rows = []

for species, names in sorted(species_to_tips.items()):

    species_set = frozenset(names)
    n = len(names)

    # Singletons cannot meaningfully test monophyly
    if n < 2:
        status = "singleton"
    else:
        status = (
            "monophyletic"
            if species_set in splits
            else "non_monophyletic"
        )

    rows.append({
        "species": species,
        "n_sequences": n,
        "status": status
    })

species_df = pd.DataFrame(rows)

testable = species_df[
    species_df["status"] != "singleton"
]

mono = (
    testable["status"] == "monophyletic"
).sum()

nonmono = (
    testable["status"] == "non_monophyletic"
).sum()

n_tested = len(testable)

mono_pct = (
    mono / n_tested * 100
    if n_tested else 0
)

# --------------------------------------------------
# Terminal branch-length QC
# --------------------------------------------------

branch_rows = []

for tip in tips:

    branch_rows.append({
        "tip": tip.name,
        "species": get_species(tip),
        "branch_length": (
            tip.branch_length
            if tip.branch_length is not None
            else 0
        )
    })

branches = pd.DataFrame(branch_rows)

q1 = branches["branch_length"].quantile(0.25)
median = branches["branch_length"].median()
q3 = branches["branch_length"].quantile(0.75)

iqr = q3 - q1

# conservative long-branch flag
threshold = q3 + 3 * iqr

branches["long_branch_outlier"] = (
    branches["branch_length"] > threshold
)

# --------------------------------------------------
# Save results
# --------------------------------------------------

species_df.to_csv(
    out_dir / f"{args.marker}_species_monophyly.csv",
    index=False
)

branches.to_csv(
    out_dir / f"{args.marker}_terminal_branches.csv",
    index=False
)

nonmono_df = species_df[
    species_df["status"] == "non_monophyletic"
]

nonmono_df.to_csv(
    out_dir / f"{args.marker}_non_monophyletic_species.csv",
    index=False
)

outliers = branches[
    branches["long_branch_outlier"]
].sort_values(
    "branch_length",
    ascending=False
)

outliers.to_csv(
    out_dir / f"{args.marker}_long_branch_outliers.csv",
    index=False
)

# --------------------------------------------------
# Summary
# --------------------------------------------------

print("\n" + "=" * 60)
print(f"{args.marker} PHYLOGENETIC SUMMARY")
print("=" * 60)

print(f"Tree tips:                {len(tips)}")
print(f"Species represented:      {len(species_df)}")
print(f"Species tested:           {n_tested}")
print(f"Monophyletic species:     {mono}")
print(f"Non-monophyletic species: {nonmono}")
print(f"Monophyly:                {mono_pct:.1f}%")

print("\nTerminal branch lengths")
print(f"Q1:                       {q1:.6f}")
print(f"Median:                   {median:.6f}")
print(f"Q3:                       {q3:.6f}")
print(f"Long-branch threshold:    {threshold:.6f}")
print(
    f"Long-branch outliers:     "
    f"{branches['long_branch_outlier'].sum()}"
)

if len(outliers):
    print("\nLongest suspicious tips:")
    print(
        outliers[
            ["tip", "branch_length"]
        ]
        .head(10)
        .to_string(index=False)
    )

print(f"\nResults saved to: {out_dir}")