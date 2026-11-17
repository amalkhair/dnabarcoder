from pathlib import Path
import pandas as pd

CLASS_DIR = Path("input/publication_scale_classification")
COVERAGE_FILE = Path(
    "results/genus_barcode_coverage/genus_barcode_coverage_long.csv"
)

OUT_DIR = Path("results/genus_barcode_overlap")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_MARKER_SEQUENCES = 200


# --------------------------------------------------
# LOAD COVERAGE AND KEEP >=200 ALTERNATIVE MARKERS
# --------------------------------------------------

coverage = pd.read_csv(COVERAGE_FILE)

candidates = coverage[
    (coverage["marker"] != "ITS")
    & (coverage["sequence_count"] >= MIN_MARKER_SEQUENCES)
].copy()


# --------------------------------------------------
# LOAD CLASSIFICATION FILES
# --------------------------------------------------

classifications = {}

for file in CLASS_DIR.glob("*.classification"):

    marker = file.stem

    df = pd.read_csv(
        file,
        sep="\t",
        dtype=str
    )

    classifications[marker] = df


its_df = classifications["ITS"]


# --------------------------------------------------
# COMPARE EACH MARKER WITH ITS
# --------------------------------------------------

rows = []

for row in candidates.itertuples():

    genus = row.genus
    marker = row.marker

    marker_df = classifications[marker]

    its_genus = its_df[
        its_df["genus"] == genus
    ].copy()

    marker_genus = marker_df[
        marker_df["genus"] == genus
    ].copy()

    # Species counts per marker
    its_counts = (
        its_genus
        .groupby("species")
        .size()
    )

    marker_counts = (
        marker_genus
        .groupby("species")
        .size()
    )

    its_species = set(its_counts.index)
    marker_species = set(marker_counts.index)

    shared_species = (
        its_species & marker_species
    )

    # Species with >=2 sequences in BOTH markers
    shared_ge2 = [
        species
        for species in shared_species
        if its_counts[species] >= 2
        and marker_counts[species] >= 2
    ]

    display_marker = (
        "CaM" if marker == "CAL" else marker
    )

    rows.append({
        "genus": genus,
        "marker": marker,
        "display_marker": display_marker,

        "ITS_sequences": len(its_genus),
        "ITS_species": len(its_species),

        "marker_sequences": len(marker_genus),
        "marker_species": len(marker_species),

        "shared_species": len(shared_species),

        "shared_species_percent_of_ITS": (
            len(shared_species)
            / len(its_species) * 100
            if its_species else 0
        ),

        "shared_species_ge2_both": len(shared_ge2),
    })


result = pd.DataFrame(rows)

result = result.sort_values(
    ["genus", "shared_species_ge2_both"],
    ascending=[True, False]
)


# --------------------------------------------------
# SAVE
# --------------------------------------------------

output = (
    OUT_DIR /
    "ITS_alternative_overlap_ge200.csv"
)

result.to_csv(
    output,
    index=False
)


# --------------------------------------------------
# PRINT
# --------------------------------------------------

print("\nITS ↔ ALTERNATIVE BARCODE OVERLAP")
print("=" * 110)

for genus in sorted(result["genus"].unique()):

    print(f"\n{genus}")
    print("-" * 110)

    subset = result[
        result["genus"] == genus
    ]

    for r in subset.itertuples():

        print(
            f"{r.display_marker:5} | "
            f"{r.marker_sequences:5} marker seq | "
            f"{r.marker_species:4} marker spp | "
            f"{r.shared_species:4} shared with ITS | "
            f"{r.shared_species_ge2_both:4} shared spp >=2 seq in both"
        )


print("\nSaved:")
print(output)