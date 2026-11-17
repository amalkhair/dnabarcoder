from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import Phylo


TREE_FILE = Path(
    "results/phylo_validation/fusarium_balanced/trees/CAL/"
    "Fusarium_CAL_matched_balanced.mafft.aligned.fas.treefile"
)

OUTPUT_DIR = Path(
    "results/phylo_validation/fusarium_balanced/tree_analysis"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------
# HELPERS
# ----------------------------------------------------------

def species_from_tip(tip_name):
    """
    FASTA/tree labels look like:

    Fusarium_species_name|CAL|rep1|ACCESSION

    Return:
    Fusarium species name
    """

    species = tip_name.split("|")[0]
    return species.replace("_", " ")


# ----------------------------------------------------------
# LOAD TREE
# ----------------------------------------------------------

tree = Phylo.read(TREE_FILE, "newick")

terminals = tree.get_terminals()

print(f"Tree tips: {len(terminals)}")


# ----------------------------------------------------------
# GROUP TIPS BY SPECIES
# ----------------------------------------------------------

species_tips = defaultdict(list)

for tip in terminals:
    species = species_from_tip(tip.name)
    species_tips[species].append(tip)


print(f"Species represented: {len(species_tips)}")


# ----------------------------------------------------------
# CREATE ALL BIPARTITIONS
#
# IQ-TREE outputs an unrooted tree. Testing monophyly only
# using common_ancestor() can depend on the arbitrary root.
#
# Instead, for every branch we record both sides of the
# corresponding split. A species is monophyletic if one
# side of any split contains exactly all of its tips.
# ----------------------------------------------------------

all_tip_names = {
    tip.name
    for tip in terminals
}

splits = set()

for clade in tree.find_clades(order="postorder"):

    descendants = {
        tip.name
        for tip in clade.get_terminals()
    }

    if (
        len(descendants) == 0
        or len(descendants) == len(all_tip_names)
    ):
        continue

    complement = all_tip_names - descendants

    splits.add(frozenset(descendants))
    splits.add(frozenset(complement))


# ----------------------------------------------------------
# TEST SPECIES MONOPHYLY
# ----------------------------------------------------------

rows = []

for species, tips in sorted(species_tips.items()):

    tip_names = {
        tip.name
        for tip in tips
    }

    n_sequences = len(tip_names)

    # A singleton cannot meaningfully test species clustering.
    if n_sequences < 2:

        rows.append({
            "species": species,
            "n_sequences": n_sequences,
            "testable": False,
            "monophyletic": None,
        })

        continue

    monophyletic = (
        frozenset(tip_names)
        in splits
    )

    rows.append({
        "species": species,
        "n_sequences": n_sequences,
        "testable": True,
        "monophyletic": monophyletic,
    })


monophyly_df = pd.DataFrame(rows)


# ----------------------------------------------------------
# SUMMARY
# ----------------------------------------------------------

testable = monophyly_df[
    monophyly_df["testable"]
].copy()

n_testable = len(testable)

n_monophyletic = int(
    testable["monophyletic"].sum()
)

n_non_monophyletic = (
    n_testable - n_monophyletic
)

percent_monophyletic = (
    n_monophyletic / n_testable * 100
    if n_testable > 0
    else np.nan
)


print("\n" + "=" * 60)
print("SPECIES MONOPHYLY")
print("=" * 60)

print(
    f"Total species: {len(monophyly_df)}"
)

print(
    f"Singleton species excluded: "
    f"{(~monophyly_df['testable']).sum()}"
)

print(
    f"Species with >=2 sequences: {n_testable}"
)

print(
    f"Monophyletic species: {n_monophyletic}"
)

print(
    f"Non-monophyletic species: {n_non_monophyletic}"
)

print(
    f"Monophyletic percentage: "
    f"{percent_monophyletic:.1f}%"
)


# ----------------------------------------------------------
# SAVE SPECIES RESULTS
# ----------------------------------------------------------

monophyly_file = (
    OUTPUT_DIR /
    "CAL_species_monophyly.csv"
)

monophyly_df.to_csv(
    monophyly_file,
    index=False
)


nonmono_file = (
    OUTPUT_DIR /
    "CAL_non_monophyletic_species.csv"
)

testable[
    testable["monophyletic"] == False
].to_csv(
    nonmono_file,
    index=False
)


# ----------------------------------------------------------
# TERMINAL BRANCH-LENGTH OUTLIERS
# ----------------------------------------------------------

branch_rows = []

for tip in terminals:

    length = (
        tip.branch_length
        if tip.branch_length is not None
        else 0
    )

    branch_rows.append({
        "tip": tip.name,
        "species": species_from_tip(tip.name),
        "terminal_branch_length": length,
    })


branch_df = pd.DataFrame(branch_rows)

q1 = branch_df[
    "terminal_branch_length"
].quantile(0.25)

q3 = branch_df[
    "terminal_branch_length"
].quantile(0.75)

iqr = q3 - q1

# Conservative threshold because real fungal data
# can naturally contain variable branch lengths.
outlier_threshold = q3 + (3 * iqr)

branch_df["long_branch_outlier"] = (
    branch_df["terminal_branch_length"]
    > outlier_threshold
)


print("\n" + "=" * 60)
print("TERMINAL BRANCH LENGTHS")
print("=" * 60)

print(
    f"Median terminal branch: "
    f"{branch_df['terminal_branch_length'].median():.4f}"
)

print(
    f"Q1: {q1:.4f}"
)

print(
    f"Q3: {q3:.4f}"
)

print(
    f"Long-branch threshold (Q3 + 3×IQR): "
    f"{outlier_threshold:.4f}"
)

print(
    f"Potential long-branch outliers: "
    f"{branch_df['long_branch_outlier'].sum()}"
)


print("\nLongest terminal branches:\n")

print(
    branch_df.sort_values(
        "terminal_branch_length",
        ascending=False
    )
    .head(15)[
        [
            "species",
            "tip",
            "terminal_branch_length",
            "long_branch_outlier"
        ]
    ]
    .to_string(index=False)
)


branch_file = (
    OUTPUT_DIR /
    "CAL_terminal_branch_lengths.csv"
)

branch_df.to_csv(
    branch_file,
    index=False
)


# ----------------------------------------------------------
# SAVE OVERALL SUMMARY
# ----------------------------------------------------------

summary = pd.DataFrame([{
    "marker": "CAL",
    "tree_tips": len(terminals),
    "total_species": len(monophyly_df),
    "singleton_species": int(
        (~monophyly_df["testable"]).sum()
    ),
    "testable_species": n_testable,
    "monophyletic_species": n_monophyletic,
    "non_monophyletic_species": n_non_monophyletic,
    "percent_monophyletic": percent_monophyletic,
    "long_branch_outliers": int(
        branch_df["long_branch_outlier"].sum()
    ),
}])


summary_file = (
    OUTPUT_DIR /
    "CAL_tree_summary.csv"
)

summary.to_csv(
    summary_file,
    index=False
)


print(f"\nSaved: {summary_file}")
print(f"Saved: {monophyly_file}")
print(f"Saved: {nonmono_file}")
print(f"Saved: {branch_file}")