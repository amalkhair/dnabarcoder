from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO


# ============================================================
# PATHS
# ============================================================

FASTA_DIR = Path("input/publication_scale_clean")
CLASS_DIR = Path("input/publication_scale_classification")

OVERLAP_FILE = Path(
    "results/genus_barcode_overlap/ITS_alternative_overlap_ge200.csv"
)

OUT_ROOT = Path("results/phyologentic_datasets")
OUT_ROOT.mkdir(parents=True, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def display_marker(marker):
    return "CaM" if marker == "CAL" else marker


def load_genus_records(marker, genus):
    """
    Return:
        species -> list of FASTA records

    Uses the classification file to connect accession IDs
    to genus/species names.
    """

    class_file = CLASS_DIR / f"{marker}.classification"
    fasta_file = FASTA_DIR / f"{marker}.fasta"

    classification = pd.read_csv(
        class_file,
        sep="\t",
        dtype=str
    )

    classification = classification[
        classification["genus"] == genus
    ].copy()

    # accession -> species
    id_to_species = dict(
        zip(
            classification["id"],
            classification["species"]
        )
    )

    records_by_species = defaultdict(list)

    fasta_ids_found = set()

    for record in SeqIO.parse(fasta_file, "fasta"):

        accession = record.id

        if accession not in id_to_species:
            continue

        species = id_to_species[accession]

        records_by_species[species].append(record)
        fasta_ids_found.add(accession)

    return records_by_species


def write_records(
    records_by_species,
    species_to_keep,
    marker_name,
    output_file
):
    """
    Write all sequences belonging to retained species.

    FASTA headers:
    Genus_species|marker|accession
    """

    output_records = []

    for species in sorted(species_to_keep):

        records = records_by_species[species]

        for record in records:

            new_record = record[:]

            safe_species = species.replace(" ", "_")

            new_record.id = (
                f"{safe_species}|"
                f"{marker_name}|"
                f"{record.id}"
            )

            new_record.name = new_record.id
            new_record.description = ""

            output_records.append(new_record)

    SeqIO.write(
        output_records,
        output_file,
        "fasta"
    )

    return len(output_records)


# ============================================================
# LOAD >=200 CANDIDATES
# ============================================================

overlap = pd.read_csv(OVERLAP_FILE)

manifest_rows = []


# ============================================================
# BUILD EACH ITS-vs-MARKER DATASET
# ============================================================

for row in overlap.itertuples():

    genus = row.genus
    marker = row.marker
    shown_marker = display_marker(marker)

    print("\n" + "=" * 70)
    print(f"{genus}: ITS vs {shown_marker}")
    print("=" * 70)

    # Load sequences using actual FASTA + classification
    its_records = load_genus_records(
        "ITS",
        genus
    )

    marker_records = load_genus_records(
        marker,
        genus
    )

    # Species must have >=2 actual sequences in BOTH datasets
    shared_species = sorted(
        species
        for species in (
            set(its_records)
            & set(marker_records)
        )
        if len(its_records[species]) >= 2
        and len(marker_records[species]) >= 2
    )

    # --------------------------------------------------------
    # OUTPUT DIRECTORY
    # --------------------------------------------------------

    pair_dir = (
        OUT_ROOT
        / genus
        / shown_marker
    )

    pair_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # SPECIES LIST
    # --------------------------------------------------------

    species_rows = []

    for species in shared_species:

        species_rows.append({
            "species": species,
            "ITS_sequences": len(
                its_records[species]
            ),
            f"{shown_marker}_sequences": len(
                marker_records[species]
            ),
        })

    species_df = pd.DataFrame(
        species_rows
    )

    species_csv = (
        pair_dir
        / f"{genus}_ITS_vs_{shown_marker}_species.csv"
    )

    species_df.to_csv(
        species_csv,
        index=False
    )

    # Simple text file as well
    species_txt = (
        pair_dir
        / f"{genus}_ITS_vs_{shown_marker}_species.txt"
    )

    with open(
        species_txt,
        "w"
    ) as handle:

        for species in shared_species:
            handle.write(
                species + "\n"
            )

    # --------------------------------------------------------
    # FASTA FILES
    # --------------------------------------------------------

    its_fasta = (
        pair_dir
        / f"{genus}_ITS_for_{shown_marker}.fasta"
    )

    marker_fasta = (
        pair_dir
        / f"{genus}_{shown_marker}_for_ITS.fasta"
    )

    n_its = write_records(
        its_records,
        shared_species,
        "ITS",
        its_fasta
    )

    n_marker = write_records(
        marker_records,
        shared_species,
        shown_marker,
        marker_fasta
    )

    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    expected_species = int(
        row.shared_species_ge2_both
    )

    status = (
        "OK"
        if len(shared_species)
        == expected_species
        else "CHECK"
    )

    print(
        f"Expected shared species: {expected_species}"
    )

    print(
        f"FASTA-confirmed species: {len(shared_species)}"
    )

    print(
        f"ITS sequences written:   {n_its}"
    )

    print(
        f"{shown_marker} sequences written: "
        f"{n_marker}"
    )

    print(
        f"QC: {status}"
    )

    manifest_rows.append({
        "genus": genus,
        "marker": shown_marker,
        "expected_shared_species_ge2": (
            expected_species
        ),
        "fasta_confirmed_shared_species": (
            len(shared_species)
        ),
        "ITS_sequences_written": n_its,
        "marker_sequences_written": n_marker,
        "QC": status,
        "directory": str(pair_dir),
    })


# ============================================================
# MASTER MANIFEST
# ============================================================

manifest = pd.DataFrame(
    manifest_rows
)

manifest = manifest.sort_values(
    [
        "genus",
        "fasta_confirmed_shared_species"
    ],
    ascending=[True, False]
)

manifest_file = (
    OUT_ROOT
    / "dataset_manifest.csv"
)

manifest.to_csv(
    manifest_file,
    index=False
)


print("\n\n" + "=" * 70)
print("DONE")
print("=" * 70)

print(
    f"\nDatasets saved to:\n{OUT_ROOT}"
)

print(
    f"\nManifest:\n{manifest_file}"
)

print("\nQC SUMMARY")
print(
    manifest[
        [
            "genus",
            "marker",
            "fasta_confirmed_shared_species",
            "ITS_sequences_written",
            "marker_sequences_written",
            "QC",
        ]
    ].to_string(index=False)
)