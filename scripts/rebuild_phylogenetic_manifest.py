from pathlib import Path
import pandas as pd
from Bio import SeqIO


ROOT = Path("results/phylogenetic_datasets")

OVERLAP_FILE = Path(
    "results/genus_barcode_overlap/ITS_alternative_overlap_ge200.csv"
)

COVERAGE_FILE = Path(
    "results/genus_barcode_coverage/genus_barcode_coverage_long.csv"
)

OUT = ROOT / "phylogenetic_dataset_manifest.csv"


def count_fasta(path):
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


coverage = pd.read_csv(COVERAGE_FILE)
overlap = pd.read_csv(OVERLAP_FILE)

coverage["display_marker"] = coverage["marker"].replace(
    {"CAL": "CaM"}
)

overlap["display_marker"] = overlap["marker"].replace(
    {"CAL": "CaM"}
)

rows = []

for row in overlap.itertuples():

    genus = row.genus
    marker = row.display_marker

    folder = ROOT / genus / marker

    species_csv = (
        folder /
        f"{genus}_ITS_vs_{marker}_species.csv"
    )

    species_txt = (
        folder /
        f"{genus}_ITS_vs_{marker}_species.txt"
    )

    its_fasta = (
        folder /
        f"{genus}_ITS_for_{marker}.fasta"
    )

    alt_fasta = (
        folder /
        f"{genus}_{marker}_for_ITS.fasta"
    )

    # Coverage information
    cov = coverage[
        (coverage["genus"] == genus)
        & (coverage["display_marker"] == marker)
    ]

    if cov.empty:
        barcode_sequences = None
        barcode_species = None
    else:
        cov = cov.iloc[0]
        barcode_sequences = int(cov["sequence_count"])
        barcode_species = int(cov["species_count"])

    # Confirm files exist
    files_ok = all(
        p.exists()
        for p in [
            species_csv,
            species_txt,
            its_fasta,
            alt_fasta,
        ]
    )

    if files_ok:
        species_df = pd.read_csv(species_csv)

        confirmed_species = len(species_df)
        its_sequences = count_fasta(its_fasta)
        alt_sequences = count_fasta(alt_fasta)

        qc = (
            "OK"
            if confirmed_species
            == int(row.shared_species_ge2_both)
            else "CHECK"
        )
    else:
        confirmed_species = None
        its_sequences = None
        alt_sequences = None
        qc = "MISSING FILE"

    rows.append({
        "genus": genus,
        "alternative_barcode": marker,

        "barcode_sequences_available":
            barcode_sequences,

        "barcode_species_available":
            barcode_species,

        "species_shared_with_ITS":
            int(row.shared_species),

        "species_shared_ge2_sequences_both":
            int(row.shared_species_ge2_both),

        "ITS_sequences_in_FASTA":
            its_sequences,

        "barcode_sequences_in_FASTA":
            alt_sequences,

        # IMPORTANT:
        # relative paths from the package root
        "species_list_csv":
            str(species_csv.relative_to(ROOT)),

        "species_list_txt":
            str(species_txt.relative_to(ROOT)),

        "ITS_FASTA":
            str(its_fasta.relative_to(ROOT)),

        "alternative_FASTA":
            str(alt_fasta.relative_to(ROOT)),

        "QC": qc,
    })


df = pd.DataFrame(rows)

df = df.sort_values(
    [
        "genus",
        "species_shared_ge2_sequences_both"
    ],
    ascending=[True, False]
)

df.to_csv(
    OUT,
    index=False
)


print("\nPHYLOGENETIC DATASET MANIFEST")
print("=" * 100)

print(
    df[
        [
            "genus",
            "alternative_barcode",
            "species_shared_ge2_sequences_both",
            "ITS_sequences_in_FASTA",
            "barcode_sequences_in_FASTA",
            "QC",
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUT)