from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = Path(
    "results/final_summary/all_marker_genus_cutoffs.csv"
)

OUTPUT_DIR = Path("results/figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Genera selected because ITS is reported in the literature
# to provide insufficient or limited species-level resolution
# in at least part of the genus.
target_genera = [
    "Aspergillus",
    "Penicillium",
    "Fusarium",
    "Trichoderma",
    "Colletotrichum",
    "Alternaria",
    "Cladosporium",
]


# Alternative biomarkers to investigate.
#
# CAL and ACT are included here even though they were excluded
# from the broad comparison, because they may be useful
# taxon-specific secondary markers.
markers = [
    "RPB1",
    "RPB2",
    "TEF1",
    "TUB2",
    "CAL",
    "ACT",
    "SSU",
    "LSU",
]


df = pd.read_csv(INPUT_FILE)

df = df[
    df["Dataset"].isin(target_genera)
    & df["marker"].isin(markers)
].copy()


# DNABarcoder stores optimized F-measure
# in the column named "confidence".
matrix = df.pivot_table(
    index="Dataset",
    columns="marker",
    values="confidence",
    aggfunc="first"
)


# Keep desired genus and biomarker order.
existing_genera = [
    genus
    for genus in target_genera
    if genus in matrix.index
]

matrix = matrix.reindex(
    index=existing_genera,
    columns=markers
)


print("\nGenus-specific optimized F-measures:\n")
print(matrix.round(3).to_string())


# ----------------------------------------------------------
# Determine best marker for every genus
# ----------------------------------------------------------

best_rows = []

for genus in matrix.index:

    available = matrix.loc[genus].dropna()

    if len(available) == 0:
        continue

    best_marker = available.idxmax()
    best_f = available.max()

    best_rows.append({
        "genus": genus,
        "best_marker": best_marker,
        "best_f_measure": best_f,
        "n_markers_available": len(available)
    })


best_df = pd.DataFrame(best_rows)

print("\nBest available marker per genus:\n")
print(best_df.round(3).to_string(index=False))


# Save underlying data
matrix.to_csv(
    OUTPUT_DIR /
    "its_problematic_genera_fmeasure_matrix.csv"
)

best_df.to_csv(
    OUTPUT_DIR /
    "its_problematic_genera_best_markers.csv",
    index=False
)


# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(11, 7)
)

masked = np.ma.masked_invalid(
    matrix.values.astype(float)
)

im = ax.imshow(
    masked,
    aspect="auto",
    vmin=0,
    vmax=1
)


ax.set_xticks(
    np.arange(len(matrix.columns))
)

ax.set_xticklabels(
    matrix.columns,
    rotation=45,
    ha="right"
)


ax.set_yticks(
    np.arange(len(matrix.index))
)

ax.set_yticklabels(
    matrix.index
)


ax.set_xlabel(
    "Alternative biomarker"
)

ax.set_ylabel(
    "Fungal genus"
)

ax.set_title(
    "Resolving performance of alternative biomarkers\n"
    "in fungal genera with reported ITS limitations"
)


# ----------------------------------------------------------
# Add values to cells
# ----------------------------------------------------------

for i, genus in enumerate(matrix.index):

    row = matrix.loc[genus]
    available = row.dropna()

    if len(available) > 0:
        best_value = available.max()

    for j, marker in enumerate(matrix.columns):

        value = matrix.iloc[i, j]

        if pd.isna(value):
            ax.text(
                j,
                i,
                "–",
                ha="center",
                va="center",
                fontsize=9
            )

        else:

            # Add * to the highest-performing available
            # marker for that genus.
            label = f"{value:.2f}"

            if value == best_value:
                label += "*"

            ax.text(
                j,
                i,
                label,
                ha="center",
                va="center",
                fontsize=9
            )


cbar = fig.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    "Optimized genus-level F-measure"
)


fig.text(
    0.5,
    0.01,
    "* Highest F-measure among available alternative markers for that genus. "
    "Dash indicates no predicted cutoff available.",
    ha="center",
    fontsize=9
)


fig.tight_layout(
    rect=[0, 0.04, 1, 1]
)


OUTPUT_FILE = (
    OUTPUT_DIR /
    "its_problematic_genera_alternative_markers.png"
)

fig.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


print(
    f"\nSaved figure: {OUTPUT_FILE}"
)