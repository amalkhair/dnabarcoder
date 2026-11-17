from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------------------------------
# INPUTS
# --------------------------------------------------

INDISTINCT_FILE = Path(
    "results/ITS_focus/dnabarcoder_cluster/ITS.indistinct.species"
)

CUTOFF_FILE = Path(
    "results/final_summary/all_marker_genus_cutoffs.csv"
)

OUT_DIR = Path("results/genus_selection")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_SEQUENCES = 1000


# --------------------------------------------------
# LOAD ITS 100%-IDENTITY GROUPS
# --------------------------------------------------

with open(INDISTINCT_FILE) as f:
    indistinct = json.load(f)


def get_genus(species_name):
    """
    Extract genus from species name.
    Also removes square brackets around names such as
    [Neocosmospora] rubicola.
    """
    return species_name.split()[0].strip("[]")


# DNABarcoder can report the same indistinguishable group
# under several species names. Convert each group to a
# canonical frozenset so it is counted only once.

unique_groups = set()

for key, members in indistinct.items():
    group = set(members)
    group.add(key)

    if len(group) >= 2:
        unique_groups.add(frozenset(group))


# --------------------------------------------------
# COUNT WITHIN-GENUS ITS AMBIGUITY
# --------------------------------------------------

genus_groups = {}
genus_species = {}

for group in unique_groups:

    species_by_genus = {}

    for species in group:
        genus = get_genus(species)

        species_by_genus.setdefault(genus, set()).add(species)

    # Only count a genus if at least TWO different species
    # from that genus occur in the same 100%-identity group.
    for genus, species_set in species_by_genus.items():

        if len(species_set) < 2:
            continue

        genus_groups.setdefault(genus, 0)
        genus_groups[genus] += 1

        genus_species.setdefault(genus, set())
        genus_species[genus].update(species_set)


ambiguity_rows = []

for genus in genus_groups:
    ambiguity_rows.append({
        "genus": genus,
        "indistinguishable_groups": genus_groups[genus],
        "ambiguous_species": len(genus_species[genus]),
    })

ambiguity = pd.DataFrame(ambiguity_rows)


# --------------------------------------------------
# LOAD ITS COVERAGE
# --------------------------------------------------

cutoffs = pd.read_csv(CUTOFF_FILE)

its = cutoffs[
    cutoffs["marker"].astype(str).str.upper() == "ITS"
].copy()

its = its[
    [
        "Dataset",
        "sequence number",
        "group number",
        "confidence",
        "cut-off",
    ]
].rename(
    columns={
        "Dataset": "genus",
        "sequence number": "ITS_sequences",
        "group number": "ITS_species",
        "confidence": "ITS_f_measure",
        "cut-off": "ITS_optimized_cutoff",
    }
)


# --------------------------------------------------
# MERGE COVERAGE + AMBIGUITY
# --------------------------------------------------

df = ambiguity.merge(
    its,
    on="genus",
    how="inner"
)

df["selected"] = df["ITS_sequences"] >= MIN_SEQUENCES

df = df.sort_values(
    "ITS_sequences",
    ascending=True
)


# Save complete table
df.to_csv(
    OUT_DIR / "ITS_genus_selection_full.csv",
    index=False
)

selected = df[df["selected"]].copy()

selected.to_csv(
    OUT_DIR / "ITS_selected_genera_final.csv",
    index=False
)


# --------------------------------------------------
# PRINT SELECTED GENERA
# --------------------------------------------------

print("\nSELECTED GENERA")
print("=" * 80)

print(
    selected[
        [
            "genus",
            "ITS_sequences",
            "ITS_species",
            "indistinguishable_groups",
            "ambiguous_species",
        ]
    ]
    .sort_values("ITS_sequences", ascending=False)
    .to_string(index=False)
)

print(
    f"\nSelected: {len(selected)} genera "
    f"with >= {MIN_SEQUENCES:,} ITS sequences."
)


# --------------------------------------------------
# FIGURE
# --------------------------------------------------

# To keep the figure readable, show genera with
# at least 400 ITS sequences.
plot_df = df[
    df["ITS_sequences"] >= 400
].copy()

plot_df = plot_df.sort_values(
    "ITS_sequences",
    ascending=True
)


fig, ax = plt.subplots(
    figsize=(10, max(6, len(plot_df) * 0.42))
)

# Selected genera are emphasized;
# excluded genera remain visually secondary.
colors = [
    "tab:blue" if selected else "lightgray"
    for selected in plot_df["selected"]
]

bars = ax.barh(
    plot_df["genus"],
    plot_df["ITS_sequences"],
    color=colors
)


# 1000-sequence threshold
ax.axvline(
    MIN_SEQUENCES,
    linestyle="--",
    linewidth=1.5,
    color="black"
)

ax.text(
    MIN_SEQUENCES + 35,
    len(plot_df) - 0.5,
    "1,000-sequence inclusion threshold",
    va="top",
    fontsize=9
)


# Add ambiguity information beside each bar
for bar, (_, row) in zip(bars, plot_df.iterrows()):

    x = bar.get_width()
    y = bar.get_y() + bar.get_height() / 2

    label = (
        f"{int(row['ITS_sequences']):,} seq | "
        f"{int(row['indistinguishable_groups'])} groups | "
        f"{int(row['ambiguous_species'])} spp."
    )

    ax.text(
        x + 30,
        y,
        label,
        va="center",
        fontsize=8
    )


ax.set_xlabel("Number of ITS sequences")
ax.set_ylabel("Fungal genus")

ax.set_title(
    "Selection of ITS-problematic fungal genera\n"
    "based on 100% sequence identity and ITS dataset coverage"
)

# Give labels enough room
ax.set_xlim(
    0,
    plot_df["ITS_sequences"].max() * 1.42
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

OUTPUT = OUT_DIR / "ITS_genus_selection.png"

plt.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved figure: {OUTPUT}")
print(
    "Saved table:",
    OUT_DIR / "ITS_genus_selection_full.csv"
)
print(
    "Saved selected genera:",
    OUT_DIR / "ITS_selected_genera_final.csv"
)