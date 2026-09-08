from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

markers = [
    "ACT", "CAL", "LSU", "PGK", "RPB1",
    "RPB2", "SSU", "TEF1", "TOP1", "TUB2"
]

input_dir = Path("input/publication_scale_clean")
output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

rows = []

for marker in markers:
    fasta_file = input_dir / f"{marker}.fasta"

    if not fasta_file.exists():
        print(f"Skipping {marker}: missing {fasta_file}")
        continue

    sequence = []

    with open(fasta_file, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    rows.append({
                        "biomarker": marker,
                        "length_bp": len("".join(sequence))
                    })
                    sequence = []
            elif line:
                sequence.append(line)

    if sequence:
        rows.append({
            "biomarker": marker,
            "length_bp": len("".join(sequence))
        })

df = pd.DataFrame(rows)

if df.empty:
    raise RuntimeError("No sequence lengths were extracted.")

df.to_csv(
    output_dir / "all_marker_sequence_lengths.csv",
    index=False
)

plot_data = [
    df.loc[df["biomarker"] == marker, "length_bp"].values
    for marker in markers
    if marker in df["biomarker"].unique()
]

plot_labels = [
    marker
    for marker in markers
    if marker in df["biomarker"].unique()
]

plt.figure(figsize=(11, 6))

plt.boxplot(
    plot_data,
    tick_labels=plot_labels,
    showfliers=False
)

plt.xlabel("Biomarker")
plt.ylabel("Sequence length (bp)")
plt.title("Sequence length distributions across fungal biomarkers")

plt.xticks(rotation=45)
plt.tight_layout()

output_file = output_dir / "sequence_length_distributions_boxplot.png"

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Extracted {len(df):,} sequence lengths.")
print(f"Saved data to: {output_dir / 'all_marker_sequence_lengths.csv'}")
print(f"Saved figure to: {output_file}")