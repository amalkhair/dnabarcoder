from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

RESULTS = Path("results")
CUTOFF_FILE = RESULTS / "final_summary/all_marker_genus_cutoffs.csv"
OUTPUT_DIR = RESULTS / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Candidate combinations emerging from Figure 5.
# We keep a few strong alternatives per genus instead of only
# the single highest F-measure.
candidates = {
    "Aspergillus": ["RPB1", "TEF1", "CAL"],
    "Penicillium": ["RPB1", "TEF1", "CAL"],
    "Fusarium": ["CAL", "RPB1", "TEF1"],
    "Trichoderma": ["TUB2", "RPB2", "CAL"],
    "Colletotrichum": ["CAL", "RPB1", "TUB2"],
    "Alternaria": ["ACT", "CAL", "RPB2"],
    "Cladosporium": ["TUB2", "TEF1", "CAL"],
}

cutoffs = pd.read_csv(CUTOFF_FILE)

rows = []

for genus, markers in candidates.items():

    for marker in markers:

        group_file = (
            RESULTS
            / f"{marker}_focus"
            / f"{marker}_indistinguishable_species_groups.csv"
        )

        if not group_file.exists():
            print(f"Missing: {group_file}")
            continue

        groups = pd.read_csv(group_file)

        # --------------------------------------------------
        # Find all species from this genus that participate
        # in ANY 100%-identity multi-species group.
        #
        # We scan every group rather than only relying on
        # representative_genus, because an identical cluster
        # may contain species from multiple genera.
        # --------------------------------------------------

        ambiguous_species = set()
        groups_containing_genus = 0

        for species_group in groups["species_group"].dropna():

            species = [
                x.strip()
                for x in species_group.split(";")
            ]

            target_species = [
                sp
                for sp in species
                if sp.startswith(genus + " ")
            ]

            if target_species:
                groups_containing_genus += 1
                ambiguous_species.update(target_species)

        # --------------------------------------------------
        # Get genus-level DNABarcoder information
        # --------------------------------------------------

        cutoff_row = cutoffs[
            (cutoffs["marker"] == marker)
            & (cutoffs["Dataset"] == genus)
        ]

        if cutoff_row.empty:
            continue

        cutoff_row = cutoff_row.iloc[0]

        # For species-rank prediction, "group number"
        # represents the number of species classes evaluated.
        total_species = int(cutoff_row["group number"])

        f_measure = float(cutoff_row["confidence"])

        n_ambiguous_species = len(ambiguous_species)

        ambiguity_rate = (
            n_ambiguous_species / total_species
            if total_species > 0
            else None
        )

        rows.append({
            "genus": genus,
            "marker": marker,
            "f_measure": f_measure,
            "total_species": total_species,
            "indistinguishable_groups": groups_containing_genus,
            "ambiguous_species": n_ambiguous_species,
            "ambiguity_rate": ambiguity_rate,
            "ambiguity_percent": ambiguity_rate * 100,
        })


df = pd.DataFrame(rows)

# Save exact underlying results
output_csv = (
    OUTPUT_DIR
    / "candidate_marker_species_ambiguity.csv"
)

df.to_csv(output_csv, index=False)

print("\nCandidate-marker species ambiguity:\n")

print(
    df[
        [
            "genus",
            "marker",
            "f_measure",
            "total_species",
            "ambiguous_species",
            "ambiguity_percent",
        ]
    ]
    .round(3)
    .to_string(index=False)
)


# ----------------------------------------------------------
# PLOT
# ----------------------------------------------------------

genus_order = list(candidates.keys())

fig, ax = plt.subplots(figsize=(12, 7))

x_positions = []
x_labels = []
values = []

position = 0

for genus in genus_order:

    genus_df = df[df["genus"] == genus]

    for _, row in genus_df.iterrows():

        x_positions.append(position)
        x_labels.append(
            f"{row['marker']}\n{genus}"
        )
        values.append(row["ambiguity_percent"])

        position += 1

    # gap between genera
    position += 0.7


bars = ax.bar(
    x_positions,
    values
)

ax.set_xticks(x_positions)
ax.set_xticklabels(
    x_labels,
    rotation=55,
    ha="right",
    fontsize=8
)

ax.set_ylabel(
    "Species involved in a 100%-identity\n"
    "multi-species group (%)"
)

ax.set_title(
    "Species-level ambiguity among candidate alternative biomarkers"
)

ax.set_ylim(
    0,
    max(values) * 1.18 if values else 1
)


# Add percentage labels
for bar, value in zip(bars, values):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f"{value:.1f}%",
        ha="center",
        va="bottom",
        fontsize=8
    )


fig.tight_layout()

output_figure = (
    OUTPUT_DIR
    / "candidate_marker_species_ambiguity.png"
)

fig.savefig(
    output_figure,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print(f"\nSaved data: {output_csv}")
print(f"Saved figure: {output_figure}")