from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# TEF1 only, because ITS is excluded
input_file = Path("input/publication_scale_clean/TEF1.fasta")
out_dir = Path("results/TEF1_focus")
out_dir.mkdir(parents=True, exist_ok=True)

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
    "biomarker": ["TEF1"],
    "n_sequences": [len(df)],
    "mean_length_bp": [round(df["length_bp"].mean(), 2)],
    "median_length_bp": [round(df["length_bp"].median(), 2)],
    "min_length_bp": [df["length_bp"].min()],
    "max_length_bp": [df["length_bp"].max()],
    "range_bp": [f'{df["length_bp"].min()}–{df["length_bp"].max()}'],
    "std_length_bp": [round(df["length_bp"].std(), 2)]
})

summary.to_csv(out_dir / "TEF1_sequence_length_summary.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["length_bp"], bins=50)
plt.xlabel("Sequence length (bp)")
plt.ylabel("Number of sequences")
plt.title("TEF1 sequence length distribution")
plt.tight_layout()
plt.savefig(out_dir / "TEF1_sequence_length_distribution.png", dpi=300)

print(summary)
print(f"\nSaved outputs to: {out_dir}")