from pathlib import Path
import pandas as pd


# --------------------------------------------------
# PATHS
# --------------------------------------------------

COVERAGE_FILE = Path(
    "results/genus_barcode_coverage/genus_barcode_coverage_long.csv"
)

OVERLAP_FILE = Path(
    "results/genus_barcode_overlap/ITS_alternative_overlap_ge200.csv"
)

PHYLO_MANIFEST = Path(
    "results/phyologentic_datasets/dataset_manifest.csv"
)

OUT_ROOT = Path("results/phyologentic_datasets")

OUT_FILE = OUT_ROOT / "handoff_manifest.csv"


# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

coverage = pd.read_csv(COVERAGE_FILE)
overlap = pd.read_csv(OVERLAP_FILE)
phylo = pd.read_csv(PHYLO_MANIFEST)


# --------------------------------------------------
# NORMALIZE MARKER NAMES
# --------------------------------------------------

coverage["display_marker"] = coverage["marker"].replace(
    {"CAL": "CaM"}
)

overlap["display_marker"] = overlap["marker"].replace(
    {"CAL": "CaM"}
)


# --------------------------------------------------
# KEEP ONLY >=200 CANDIDATES
# --------------------------------------------------

coverage = coverage[
    (coverage["marker"] != "ITS")
    & (coverage["sequence_count"] >= 200)
].copy()


# --------------------------------------------------
# MERGE COVERAGE + OVERLAP
# --------------------------------------------------

summary = overlap.merge(
    coverage[
        [
            "genus",
            "display_marker",
            "sequence_count",
            "species_count",
        ]
    ],
    on=["genus", "display_marker"],
    how="left",
)


# --------------------------------------------------
# MERGE FASTA MANIFEST
# --------------------------------------------------

summary = summary.merge(
    phylo[
        [
            "genus",
            "marker",
            "fasta_confirmed_shared_species",
            "ITS_sequences_written",
            "marker_sequences_written",
            "QC",
        ]
    ],
    left_on=["genus", "display_marker"],
    right_on=["genus", "marker"],
    how="left",
)


# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

def make_paths(row):

    genus = row["genus"]
    marker = row["display_marker"]

    base = Path(genus) / marker

    return pd.Series({
        "species_list_csv": str(
            base /
            f"{genus}_ITS_vs_{marker}_species.csv"
        ),

        "species_list_txt": str(
            base /
            f"{genus}_ITS_vs_{marker}_species.txt"
        ),

        "ITS_FASTA": str(
            base /
            f"{genus}_ITS_for_{marker}.fasta"
        ),

        "alternative_FASTA": str(
            base /
            f"{genus}_{marker}_for_ITS.fasta"
        ),
    })


paths = summary.apply(
    make_paths,
    axis=1
)

summary = pd.concat(
    [summary, paths],
    axis=1
)


# --------------------------------------------------
# CLEAN FINAL TABLE
# --------------------------------------------------

final = pd.DataFrame({
    "genus": summary["genus"],

    "alternative_barcode": summary[
        "display_marker"
    ],

    "barcode_sequences_available": summary[
        "sequence_count"
    ],

    "barcode_species_available": summary[
        "species_count"
    ],

    "species_shared_with_ITS": summary[
        "shared_species"
    ],

    "species_shared_ge2_sequences_both": summary[
        "shared_species_ge2_both"
    ],

    "ITS_sequences_in_handoff_FASTA": summary[
        "ITS_sequences_written"
    ],

    "barcode_sequences_in_handoff_FASTA": summary[
        "marker_sequences_written"
    ],

    "species_list_csv": summary[
        "species_list_csv"
    ],

    "species_list_txt": summary[
        "species_list_txt"
    ],

    "ITS_FASTA": summary[
        "ITS_FASTA"
    ],

    "alternative_FASTA": summary[
        "alternative_FASTA"
    ],

    "QC": summary["QC"],
})


# Sort within genus by usable shared species
final = final.sort_values(
    [
        "genus",
        "species_shared_ge2_sequences_both"
    ],
    ascending=[True, False]
)


# --------------------------------------------------
# SAVE
# --------------------------------------------------

final.to_csv(
    OUT_FILE,
    index=False
)


print("\nDUONG PHYLOGENY HANDOFF")
print("=" * 100)

print(
    final[
        [
            "genus",
            "alternative_barcode",
            "barcode_sequences_available",
            "barcode_species_available",
            "species_shared_ge2_sequences_both",
            "QC",
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUT_FILE)