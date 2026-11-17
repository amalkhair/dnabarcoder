from pathlib import Path
import pandas as pd

CLASS_DIR = Path("input/publication_scale_classification")
OUT_DIR = Path("results/phylo_validation")
OUT_DIR.mkdir(parents=True, exist_ok=True)

genera = [
    "Aspergillus",
    "Penicillium",
    "Fusarium",
    "Trichoderma",
    "Colletotrichum",
    "Alternaria",
    "Cladosporium",
]

markers = [
    "RPB1",
    "RPB2",
    "TEF1",
    "TUB2",
    "CAL",   # CaM
    "ACT",
]

# --------------------------------------------------
# Load all classification files once
# --------------------------------------------------

data = {}

for marker in ["ITS"] + markers:
    file = CLASS_DIR / f"{marker}.classification"

    if not file.exists():
        print(f"Missing: {file}")
        continue

    df = pd.read_csv(file, sep="\t")

    df = df.dropna(subset=["genus", "species"]).copy()
    df["genus"] = df["genus"].astype(str).str.strip()
    df["species"] = df["species"].astype(str).str.strip()

    data[marker] = df


# --------------------------------------------------
# Compare ITS with every alternative marker
# --------------------------------------------------

rows = []

for genus in genera:

    its = data["ITS"]
    its_g = its[its["genus"] == genus]

    its_species = set(its_g["species"].unique())

    for marker in markers:

        if marker not in data:
            continue

        alt = data[marker]
        alt_g = alt[alt["genus"] == genus]

        alt_species = set(alt_g["species"].unique())

        shared = its_species & alt_species

        # Sequence counts for shared species only
        its_shared = its_g[
            its_g["species"].isin(shared)
        ]

        alt_shared = alt_g[
            alt_g["species"].isin(shared)
        ]

        its_counts = its_shared.groupby("species").size()
        alt_counts = alt_shared.groupby("species").size()

        # Species where BOTH markers have >=2 sequences
        testable_species = [
            sp for sp in shared
            if its_counts.get(sp, 0) >= 2
            and alt_counts.get(sp, 0) >= 2
        ]

        rows.append({
            "genus": genus,
            "marker": marker,

            "ITS_species": len(its_species),
            "marker_species": len(alt_species),
            "shared_species": len(shared),

            "shared_percent_of_ITS": (
                len(shared) / len(its_species) * 100
                if its_species else 0
            ),

            "ITS_matched_sequences": len(its_shared),
            "marker_matched_sequences": len(alt_shared),

            "shared_species_with_2plus_in_both":
                len(testable_species),

            "ITS_median_seq_per_shared_species": (
                its_counts.median()
                if len(its_counts)
                else 0
            ),

            "marker_median_seq_per_shared_species": (
                alt_counts.median()
                if len(alt_counts)
                else 0
            ),
        })


result = pd.DataFrame(rows)

# Sort best overlaps first within each genus
result = result.sort_values(
    ["genus", "shared_species"],
    ascending=[True, False]
)

output = OUT_DIR / "ITS_candidate_species_overlap.csv"
result.to_csv(output, index=False)


# --------------------------------------------------
# Print compact result
# --------------------------------------------------

print("\nITS ↔ ALTERNATIVE MARKER OVERLAP\n")

print(
    result[
        [
            "genus",
            "marker",
            "ITS_species",
            "marker_species",
            "shared_species",
            "shared_percent_of_ITS",
            "shared_species_with_2plus_in_both",
        ]
    ]
    .round(1)
    .to_string(index=False)
)

print(f"\nSaved: {output}")