from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path(
    "results/final_summary/all_marker_genus_cutoffs.csv"
)

output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

markers = [
    "RPB1",
    "RPB2",
    "TEF1",
    "TUB2",
    "SSU",
    "LSU"
]

df = pd.read_csv(input_file)

df = df[df["marker"].isin(markers)].copy()

# DNABarcoder's "confidence" field contains
# the optimized F-measure.
matrix = df.pivot_table(
    index="Dataset",
    columns="marker",
    values="confidence",
    aggfunc="first"
)

# Keep only genera that have results for ALL six markers.
matrix = matrix.dropna(subset=markers)

# Keep desired marker order.
matrix = matrix[markers]

# Alphabetical genus order for transparency.
matrix = matrix.sort_index()

print(f"Genera represented by all six markers: {len(matrix)}")
print(matrix.round(3).to_string())

# Save exact underlying data.
matrix.to_csv(
    output_dir / "shared_genera_fmeasure_matrix.csv"
)

fig, ax = plt.subplots(figsize=(9, 8))

im = ax.imshow(
    matrix.values,
    aspect="auto",
    vmin=0,
    vmax=1
)

ax.set_xticks(np.arange(len(markers)))
ax.set_xticklabels(markers)

ax.set_yticks(np.arange(len(matrix.index)))
ax.set_yticklabels(matrix.index)

ax.set_xlabel("Biomarker")
ax.set_ylabel("Fungal genus")

ax.set_title(
    "Genus-specific resolving performance across\n"
    "well-represented alternative fungal biomarkers"
)

# Write F-measure inside every cell.
for i in range(len(matrix.index)):
    for j in range(len(markers)):
        value = matrix.iloc[i, j]

        ax.text(
            j,
            i,
            f"{value:.2f}",
            ha="center",
            va="center",
            fontsize=8
        )

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Optimized F-measure")

fig.tight_layout()

output_file = (
    output_dir /
    "shared_genera_fmeasure_heatmap.png"
)

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print(f"\nSaved: {output_file}")