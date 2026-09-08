from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage: python scripts/marker_length_stats.py MARKER")
    print("Example: python scripts/marker_length_stats.py TEF1")
    sys.exit(1)

marker = sys.argv[1].upper()

input_file = Path(f"input/publication_scale_clean/{marker}.fasta")
out_dir = Path(f"results/{marker}_focus")
out_dir.mkdir(parents=True, exist_ok=True)

if not input_file.exists():
    raise FileNotFoundError(f"Missing input FASTA: {input_file}")

lengths = []
seq = ""

with open(input_file, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if seq:
                lengths.append(len(seq))
                seq = ""
        else:
            seq += line

if seq:
    lengths.append(len(seq))

df = pd.DataFrame({"length_bp": lengths})

summary = pd.DataFrame({
    "biomarker": [marker],
    "n_sequences": [len(df)],
    "mean_length_bp": [round(df["length_bp"].mean(), 2)],
    "median_length_bp": [round(df["length_bp"].median(), 2)],
    "min_length_bp": [df["length_bp"].min()],
    "max_length_bp": [df["length_bp"].max()],
    "range_bp": [f'{df["length_bp"].min()}–{df["length_bp"].max()}'],
    "std_length_bp": [round(df["length_bp"].std(), 2)]
})

summary.to_csv(out_dir / f"{marker}_sequence_length_summary.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["length_bp"], bins=50)
plt.xlabel("Sequence length (bp)")
plt.ylabel("Number of sequences")
plt.title(f"{marker} sequence length distribution")
plt.tight_layout()
plt.savefig(out_dir / f"{marker}_sequence_length_distribution.png", dpi=300)

print(summary)
print(f"\nSaved outputs to: {out_dir}")