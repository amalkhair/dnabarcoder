from pathlib import Path
import random
import pandas as pd
from Bio import SeqIO

GENUS = "Fusarium"
MARKERS = ["ITS", "CAL"]   # CAL = CaM
MAX_PER_SPECIES = 3
SEED = 42

FASTA_DIR = Path("input/publication_scale_clean")
CLASS_DIR = Path("input/publication_scale_classification")

OUT_DIR = Path(
    "results/phylo_validation/fusarium_ITS_vs_CaM"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

random.seed(SEED)

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

data = {}

for marker in MARKERS:

    classification = pd.read_csv(
        CLASS_DIR / f"{marker}.classification",
        sep="\t"
    )

    classification = classification[
        classification["genus"] == GENUS
    ].copy()

    classification["id"] = (
        classification["id"].astype(str).str.strip()
    )

    classification["species"] = (
        classification["species"].astype(str).str.strip()
    )

    fasta = {
        record.id: record
        for record in SeqIO.parse(
            FASTA_DIR / f"{marker}.fasta",
            "fasta"
        )
    }

    data[marker] = {
        "classification": classification,
        "fasta": fasta
    }


# --------------------------------------------------
# FIND SHARED SPECIES WITH >=2 SEQUENCES IN BOTH
# --------------------------------------------------

species_sets = []

for marker in MARKERS:
    df = data[marker]["classification"]

    counts = df.groupby("species").size()

    species_sets.append(
        set(counts[counts >= 2].index)
    )

shared_species = set.intersection(*species_sets)

print(
    f"Fusarium species with >=2 sequences "
    f"in both ITS and CaM: {len(shared_species)}"
)


# --------------------------------------------------
# BALANCED SAMPLING
# --------------------------------------------------

summary_rows = []
selected = {
    marker: []
    for marker in MARKERS
}

for species in sorted(shared_species):

    marker_ids = {}

    for marker in MARKERS:

        df = data[marker]["classification"]
        fasta = data[marker]["fasta"]

        ids = (
            df[df["species"] == species]["id"]
            .tolist()
        )

        # Keep only IDs actually present in FASTA
        ids = [
            seq_id
            for seq_id in ids
            if seq_id in fasta
        ]

        marker_ids[marker] = list(dict.fromkeys(ids))

    # Equal number in ITS and CaM
    n = min(
        MAX_PER_SPECIES,
        len(marker_ids["ITS"]),
        len(marker_ids["CAL"])
    )

    # We require >=2 so both trees can test monophyly
    if n < 2:
        continue

    for marker in MARKERS:

        chosen = random.sample(
            marker_ids[marker],
            n
        )

        for rep, seq_id in enumerate(chosen, start=1):

            record = data[marker]["fasta"][seq_id]

            safe_species = species.replace(" ", "_")

            new_id = (
                f"{safe_species}"
                f"|{marker}"
                f"|rep{rep}"
                f"|{seq_id}"
            )

            selected[marker].append(
                (new_id, str(record.seq))
            )

    summary_rows.append({
        "species": species,
        "sequences_selected_per_marker": n,
        "ITS_available": len(marker_ids["ITS"]),
        "CaM_available": len(marker_ids["CAL"])
    })


# --------------------------------------------------
# WRITE FASTA FILES
# --------------------------------------------------

for marker in MARKERS:

    display = "CaM" if marker == "CAL" else marker

    outfile = (
        OUT_DIR /
        f"Fusarium_{display}_matched_balanced.fasta"
    )

    with open(outfile, "w") as f:

        for header, sequence in selected[marker]:

            f.write(f">{header}\n")

            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

    print(
        f"{display}: {len(selected[marker])} sequences "
        f"saved to {outfile}"
    )


# --------------------------------------------------
# SAVE SUMMARY
# --------------------------------------------------

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT_DIR /
    "Fusarium_ITS_vs_CaM_sampling_summary.csv",
    index=False
)

with open(
    OUT_DIR / "Fusarium_species_used.txt",
    "w"
) as f:

    for species in summary["species"]:
        f.write(species + "\n")


print("\n" + "=" * 60)
print("FINAL MATCHED DATASET")
print("=" * 60)

print(
    f"Species retained: {len(summary)}"
)

print(
    "\nSequences selected per species:"
)

print(
    summary["sequences_selected_per_marker"]
    .value_counts()
    .sort_index()
    .to_string()
)

print(
    f"\nTotal ITS sequences: {len(selected['ITS'])}"
)

print(
    f"Total CaM sequences: {len(selected['CAL'])}"
)