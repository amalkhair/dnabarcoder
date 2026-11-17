from pathlib import Path
import random
import pandas as pd


# ----------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------

MARKERS = ["CAL", "RPB1", "TEF1"]
GENUS = "Fusarium"

MAX_SEQS_PER_SPECIES = 3
RANDOM_SEED = 42

CLASS_DIR = Path("input/publication_scale_classification")
FASTA_DIR = Path("input/publication_scale_clean")

OUT_DIR = Path("results/phylo_validation/fusarium_balanced")
OUT_DIR.mkdir(parents=True, exist_ok=True)

random.seed(RANDOM_SEED)


# ----------------------------------------------------------
# SIMPLE FASTA READER / WRITER
# ----------------------------------------------------------

def read_fasta(path):
    """
    Returns:
        dict:
        sequence_id -> (full_header, sequence)
    """

    records = {}

    header = None
    seq_parts = []

    with open(path) as f:
        for line in f:

            line = line.rstrip()

            if line.startswith(">"):

                if header is not None:
                    seq = "".join(seq_parts)

                    # Use first word of FASTA header as sequence ID
                    seq_id = header.split()[0]

                    records[seq_id] = (
                        header,
                        seq
                    )

                header = line[1:]
                seq_parts = []

            else:
                seq_parts.append(line.strip())

        # Final sequence
        if header is not None:

            seq = "".join(seq_parts)
            seq_id = header.split()[0]

            records[seq_id] = (
                header,
                seq
            )

    return records


def write_fasta(records, output_file):

    with open(output_file, "w") as f:

        for new_header, sequence in records:

            f.write(f">{new_header}\n")

            # Wrap sequence every 80 characters
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")


# ----------------------------------------------------------
# DETECT SEQUENCE-ID COLUMN
# ----------------------------------------------------------

def detect_id_column(df):

    candidates = [
        "sequence_id",
        "sequenceid",
        "seq_id",
        "seqid",
        "id",
        "ID",
        "accession",
        "accession_id",
        "identifier",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        "Could not automatically identify the sequence-ID column.\n"
        f"Columns available are:\n{list(df.columns)}"
    )


# ----------------------------------------------------------
# LOAD CLASSIFICATION + FASTA
# ----------------------------------------------------------

classification = {}
fasta_records = {}

for marker in MARKERS:

    class_file = (
        CLASS_DIR /
        f"{marker}.classification"
    )

    fasta_file = (
        FASTA_DIR /
        f"{marker}.fasta"
    )

    if not class_file.exists():
        raise FileNotFoundError(class_file)

    if not fasta_file.exists():
        raise FileNotFoundError(fasta_file)

    df = pd.read_csv(
        class_file,
        sep="\t"
    )

    id_col = detect_id_column(df)

    # Keep Fusarium only
    df = df[
        df["genus"].astype(str) == GENUS
    ].copy()

    df = df.dropna(
        subset=["species", id_col]
    )

    df["species"] = (
        df["species"]
        .astype(str)
        .str.strip()
    )

    df[id_col] = (
        df[id_col]
        .astype(str)
        .str.strip()
    )

    classification[marker] = {
        "data": df,
        "id_col": id_col,
    }

    fasta_records[marker] = read_fasta(
        fasta_file
    )

    print(
        f"{marker}: "
        f"{len(df)} Fusarium classification rows | "
        f"{df['species'].nunique()} species | "
        f"{len(fasta_records[marker])} FASTA records"
    )


# ----------------------------------------------------------
# FIND SPECIES SHARED BY ALL 3 MARKERS
# ----------------------------------------------------------

species_sets = [
    set(
        classification[m]["data"]["species"]
    )
    for m in MARKERS
]

shared_species = set.intersection(
    *species_sets
)

print(
    f"\nSpecies shared by CAL, RPB1 and TEF1: "
    f"{len(shared_species)}"
)


# ----------------------------------------------------------
# MATCH CLASSIFICATION IDs TO FASTA IDs
# ----------------------------------------------------------

species_to_ids = {}

for marker in MARKERS:

    df = classification[marker]["data"]
    id_col = classification[marker]["id_col"]

    marker_species = {}

    for species in sorted(shared_species):

        rows = df[
            df["species"] == species
        ]

        ids = []

        for seq_id in rows[id_col]:

            # Try exact ID first
            if seq_id in fasta_records[marker]:
                ids.append(seq_id)
                continue

            # Sometimes IDs contain extra whitespace/header content
            short_id = seq_id.split()[0]

            if short_id in fasta_records[marker]:
                ids.append(short_id)

        # Remove duplicates while preserving order
        ids = list(dict.fromkeys(ids))

        marker_species[species] = ids

    species_to_ids[marker] = marker_species


# ----------------------------------------------------------
# TRUE BALANCED SAMPLING
# ----------------------------------------------------------

selected = {
    marker: {}
    for marker in MARKERS
}

summary_rows = []

usable_species = []

for species in sorted(shared_species):

    available = {
        marker: len(
            species_to_ids[marker][species]
        )
        for marker in MARKERS
    }

    # Number to use must be identical in all markers
    n_to_sample = min(
        MAX_SEQS_PER_SPECIES,
        *available.values()
    )

    if n_to_sample < 1:
        print(
            f"Skipping {species}: "
            f"FASTA-linked counts = {available}"
        )
        continue

    usable_species.append(species)

    for marker in MARKERS:

        candidates = (
            species_to_ids[marker][species]
        )

        chosen = random.sample(
            candidates,
            n_to_sample
        )

        selected[marker][species] = chosen

    summary_rows.append({
        "species": species,

        "CAL_available": available["CAL"],
        "RPB1_available": available["RPB1"],
        "TEF1_available": available["TEF1"],

        "sequences_selected_per_marker":
            n_to_sample,
    })


# ----------------------------------------------------------
# WRITE BALANCED FASTAS
# ----------------------------------------------------------

for marker in MARKERS:

    output_records = []

    for species in usable_species:

        ids = selected[marker][species]

        for replicate_number, seq_id in enumerate(
            ids,
            start=1
        ):

            original_header, sequence = (
                fasta_records[marker][seq_id]
            )

            # Tree-friendly label
            safe_species = (
                species.replace(" ", "_")
            )

            new_header = (
                f"{safe_species}"
                f"|{marker}"
                f"|rep{replicate_number}"
                f"|{seq_id}"
            )

            output_records.append(
                (
                    new_header,
                    sequence
                )
            )

    output_file = (
        OUT_DIR /
        f"Fusarium_{marker}_matched_balanced.fasta"
    )

    write_fasta(
        output_records,
        output_file
    )

    print(
        f"{marker}: wrote "
        f"{len(output_records)} sequences -> "
        f"{output_file}"
    )


# ----------------------------------------------------------
# SAVE SAMPLING SUMMARY
# ----------------------------------------------------------

summary_df = pd.DataFrame(
    summary_rows
)

summary_file = (
    OUT_DIR /
    "Fusarium_balanced_sampling_summary.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)


# Save exact shared species used
species_file = (
    OUT_DIR /
    "Fusarium_species_used.txt"
)

with open(species_file, "w") as f:

    for species in usable_species:
        f.write(species + "\n")


# ----------------------------------------------------------
# FINAL CHECKS
# ----------------------------------------------------------

print("\n" + "=" * 60)
print("BALANCED DATASET SUMMARY")
print("=" * 60)

print(
    f"Shared species originally: "
    f"{len(shared_species)}"
)

print(
    f"Species successfully matched to FASTA "
    f"in all markers: {len(usable_species)}"
)

print(
    "\nSequences selected per species:"
)

print(
    summary_df[
        "sequences_selected_per_marker"
    ]
    .value_counts()
    .sort_index()
    .to_string()
)

expected_per_marker = (
    summary_df[
        "sequences_selected_per_marker"
    ].sum()
)

print(
    f"\nTotal sequences per marker: "
    f"{expected_per_marker}"
)

print(
    f"\nSaved summary: {summary_file}"
)

print(
    f"Saved species list: {species_file}"
)