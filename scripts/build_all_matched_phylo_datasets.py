from pathlib import Path
import random
import pandas as pd
from Bio import SeqIO

FASTA_DIR = Path("input/publication_scale_clean")
CLASS_DIR = Path("input/publication_scale_classification")
OUT_ROOT = Path("results/phylo_validation/final_matched_datasets")

MAX_PER_SPECIES = 3
MIN_PER_SPECIES = 2
SEED = 42

random.seed(SEED)
OUT_ROOT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# Comparisons to prepare
# CAL = CaM
# --------------------------------------------------

CASES = [
    ("Aspergillus", "CAL"),
    ("Aspergillus", "RPB1"),      # exploratory: low multi-seq overlap

    ("Penicillium", "RPB1"),
    ("Penicillium", "RPB2"),

    ("Fusarium", "CAL"),
    ("Fusarium", "RPB1"),
    ("Fusarium", "TEF1"),

    ("Trichoderma", "RPB2"),

    ("Colletotrichum", "CAL"),

    ("Alternaria", "ACT"),
    ("Alternaria", "CAL"),
    ("Alternaria", "RPB2"),

    ("Cladosporium", "TEF1"),
]


# --------------------------------------------------
# Cache datasets so files aren't repeatedly loaded
# --------------------------------------------------

cache = {}

def load_marker(marker):
    if marker in cache:
        return cache[marker]

    fasta_path = FASTA_DIR / f"{marker}.fasta"
    class_path = CLASS_DIR / f"{marker}.classification"

    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing FASTA: {fasta_path}")

    if not class_path.exists():
        raise FileNotFoundError(f"Missing classification: {class_path}")

    df = pd.read_csv(class_path, sep="\t")

    df = df.dropna(subset=["id", "genus", "species"]).copy()

    for col in ["id", "genus", "species"]:
        df[col] = df[col].astype(str).str.strip()

    fasta = {
        record.id: record
        for record in SeqIO.parse(fasta_path, "fasta")
    }

    # Only keep classification rows that actually have a sequence
    df = df[df["id"].isin(fasta)].copy()

    cache[marker] = {
        "df": df,
        "fasta": fasta
    }

    return cache[marker]


def display_marker(marker):
    return "CaM" if marker == "CAL" else marker


manifest_rows = []


# --------------------------------------------------
# Build each matched ITS-vs-marker dataset
# --------------------------------------------------

for genus, alt_marker in CASES:

    print("\n" + "=" * 70)
    print(f"{genus}: ITS vs {display_marker(alt_marker)}")
    print("=" * 70)

    its_data = load_marker("ITS")
    alt_data = load_marker(alt_marker)

    its_df = its_data["df"][
        its_data["df"]["genus"] == genus
    ].copy()

    alt_df = alt_data["df"][
        alt_data["df"]["genus"] == genus
    ].copy()

    # Sequence counts per species
    its_counts = its_df.groupby("species")["id"].nunique()
    alt_counts = alt_df.groupby("species")["id"].nunique()

    # Only species with >=2 sequences in BOTH markers
    its_testable = set(
        its_counts[its_counts >= MIN_PER_SPECIES].index
    )

    alt_testable = set(
        alt_counts[alt_counts >= MIN_PER_SPECIES].index
    )

    shared_species = sorted(
        its_testable & alt_testable
    )

    folder_name = (
        f"{genus}_ITS_vs_{display_marker(alt_marker)}"
    )

    out_dir = OUT_ROOT / folder_name
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = {
        "ITS": [],
        alt_marker: []
    }

    sampling_rows = []

    for species in shared_species:

        its_ids = (
            its_df[its_df["species"] == species]["id"]
            .drop_duplicates()
            .tolist()
        )

        alt_ids = (
            alt_df[alt_df["species"] == species]["id"]
            .drop_duplicates()
            .tolist()
        )

        # Same number of sequences in both markers
        n = min(
            MAX_PER_SPECIES,
            len(its_ids),
            len(alt_ids)
        )

        if n < MIN_PER_SPECIES:
            continue

        # Deterministic sampling
        species_seed = (
            SEED
            + sum(ord(c) for c in genus)
            + sum(ord(c) for c in species)
            + sum(ord(c) for c in alt_marker)
        )

        rng = random.Random(species_seed)

        chosen_its = rng.sample(its_ids, n)
        chosen_alt = rng.sample(alt_ids, n)

        for rep, seq_id in enumerate(chosen_its, start=1):
            record = its_data["fasta"][seq_id]

            header = (
                f"{species.replace(' ', '_')}"
                f"|ITS|rep{rep}|{seq_id}"
            )

            selected["ITS"].append(
                (header, str(record.seq))
            )

        for rep, seq_id in enumerate(chosen_alt, start=1):
            record = alt_data["fasta"][seq_id]

            header = (
                f"{species.replace(' ', '_')}"
                f"|{display_marker(alt_marker)}"
                f"|rep{rep}|{seq_id}"
            )

            selected[alt_marker].append(
                (header, str(record.seq))
            )

        sampling_rows.append({
            "genus": genus,
            "species": species,
            "alternative_marker": display_marker(alt_marker),
            "sequences_selected_per_marker": n,
            "ITS_available": len(its_ids),
            f"{display_marker(alt_marker)}_available": len(alt_ids),
        })

    # --------------------------------------------------
    # Write FASTAs
    # --------------------------------------------------

    its_out = (
        out_dir /
        f"{genus}_ITS_matched_balanced.fasta"
    )

    alt_out = (
        out_dir /
        f"{genus}_{display_marker(alt_marker)}_matched_balanced.fasta"
    )

    for marker, outfile in [
        ("ITS", its_out),
        (alt_marker, alt_out)
    ]:

        with open(outfile, "w") as f:

            for header, sequence in selected[marker]:

                f.write(f">{header}\n")

                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")

    # --------------------------------------------------
    # Save species list + sampling summary
    # --------------------------------------------------

    sampling_df = pd.DataFrame(sampling_rows)

    sampling_df.to_csv(
        out_dir / "sampling_summary.csv",
        index=False
    )

    with open(out_dir / "species_used.txt", "w") as f:
        for species in sampling_df["species"]:
            f.write(species + "\n")

    n_species = len(sampling_df)
    n_its = len(selected["ITS"])
    n_alt = len(selected[alt_marker])

    print(f"Species retained: {n_species}")
    print(f"ITS sequences:     {n_its}")
    print(
        f"{display_marker(alt_marker)} sequences: "
        f"{n_alt}"
    )

    if n_its != n_alt:
        print("WARNING: sequence counts are NOT balanced!")

    manifest_rows.append({
        "genus": genus,
        "alternative_marker": display_marker(alt_marker),
        "species_retained": n_species,
        "ITS_sequences": n_its,
        "alternative_sequences": n_alt,
        "folder": str(out_dir),
    })


# --------------------------------------------------
# Save master manifest
# --------------------------------------------------

manifest = pd.DataFrame(manifest_rows)

manifest.to_csv(
    OUT_ROOT / "dataset_manifest.csv",
    index=False
)

print("\n\n" + "=" * 70)
print("FINAL DATASET MANIFEST")
print("=" * 70)

print(
    manifest[
        [
            "genus",
            "alternative_marker",
            "species_retained",
            "ITS_sequences",
            "alternative_sequences"
        ]
    ].to_string(index=False)
)

print(
    "\nSaved all datasets to:\n"
    f"{OUT_ROOT}"
)