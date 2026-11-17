from pathlib import Path
import pandas as pd

CLASS_DIR = Path("input/publication_scale_classification")
OUT_DIR = Path("results/phylo_validation")
OUT_DIR.mkdir(parents=True, exist_ok=True)

cases = {
    "Aspergillus": ["RPB1", "TEF1"],
    "Penicillium": ["RPB1", "TEF1"],
    "Fusarium": ["CAL", "RPB1", "TEF1"],
}

summary_rows = []

for genus, markers in cases.items():

    print(f"\n{'=' * 60}")
    print(genus)
    print(f"{'=' * 60}")

    marker_data = {}
    species_sets = []

    for marker in markers:

        file = CLASS_DIR / f"{marker}.classification"

        if not file.exists():
            print(f"Missing: {file}")
            continue

        df = pd.read_csv(file, sep="\t")

        genus_df = df[df["genus"] == genus].copy()

        species = set(genus_df["species"].dropna().unique())

        marker_data[marker] = genus_df
        species_sets.append(species)

        print(
            f"{marker}: "
            f"{len(genus_df)} sequences, "
            f"{len(species)} species"
        )

    if len(species_sets) != len(markers):
        print("Skipping shared-species calculation because a marker is missing.")
        continue

    # Species represented by ALL markers in this case
    shared_species = set.intersection(*species_sets)

    print(f"\nShared species across all markers: {len(shared_species)}")

    # Save exact species list
    species_file = (
        OUT_DIR /
        f"{genus}_shared_species.txt"
    )

    with open(species_file, "w") as f:
        for species in sorted(shared_species):
            f.write(species + "\n")

    row = {
        "genus": genus,
        "markers": "+".join(markers),
        "shared_species": len(shared_species),
    }

    # Count sequences belonging specifically to the shared species
    for marker in markers:

        df = marker_data[marker]

        matched = df[
            df["species"].isin(shared_species)
        ]

        counts = matched.groupby("species").size()

        row[f"{marker}_all_species"] = (
            marker_data[marker]["species"].nunique()
        )

        row[f"{marker}_matched_sequences"] = len(matched)

        row[f"{marker}_median_sequences_per_species"] = (
            counts.median() if len(counts) else 0
        )

        row[f"{marker}_min_sequences_per_species"] = (
            counts.min() if len(counts) else 0
        )

        print(
            f"{marker} within shared set: "
            f"{len(matched)} sequences | "
            f"median {counts.median() if len(counts) else 0:.1f} "
            f"seq/species"
        )

    summary_rows.append(row)


summary = pd.DataFrame(summary_rows)

output = (
    OUT_DIR /
    "matched_species_overlap_summary.csv"
)

summary.to_csv(output, index=False)

print("\n\nSUMMARY")
print(summary.to_string(index=False))

print(f"\nSaved: {output}")