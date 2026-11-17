from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT = Path("results/final_summary/all_marker_genus_cutoffs.csv")
OUT = Path("results/figures")
OUT.mkdir(parents=True, exist_ok=True)

# Genera selected because ITS has reported species-level limitations
genera = [
    "Aspergillus",
    "Penicillium",
    "Fusarium",
    "Trichoderma",
    "Colletotrichum",
    "Alternaria",
    "Cladosporium",
]

# ITS first so comparison is immediately visible
markers = [
    "ITS",
    "RPB1",
    "RPB2",
    "TEF1",
    "TUB2",
    "CAL",
    "ACT",
    "SSU",
    "LSU",
]

display_names = {
    "ITS": "ITS",
    "RPB1": "RPB1",
    "RPB2": "RPB2",
    "TEF1": "TEF1",
    "TUB2": "TUB2",
    "CAL": "CaM/CAL",
    "ACT": "ACT",
    "SSU": "SSU",
    "LSU": "LSU",
}

df = pd.read_csv(INPUT)

sub = df[
    df["Dataset"].isin(genera) &
    df["marker"].isin(markers)
].copy()

# F-measure matrix
f_matrix = sub.pivot_table(
    index="Dataset",
    columns="marker",
    values="confidence",
    aggfunc="first"
).reindex(index=genera, columns=markers)

# Number of species evaluated
n_matrix = sub.pivot_table(
    index="Dataset",
    columns="marker",
    values="group number",
    aggfunc="first"
).reindex(index=genera, columns=markers)

# Save underlying data in long format
output_rows = []

for genus in genera:
    its_f = f_matrix.loc[genus, "ITS"]

    for marker in markers:
        f = f_matrix.loc[genus, marker]
        n = n_matrix.loc[genus, marker]

        output_rows.append({
            "genus": genus,
            "marker": marker,
            "f_measure": f,
            "n_species": n,
            "difference_vs_ITS":
                f - its_f
                if pd.notna(f) and pd.notna(its_f)
                else np.nan
        })

pd.DataFrame(output_rows).to_csv(
    OUT / "main_genus_biomarker_heatmap_data.csv",
    index=False
)

# --------------------------------------------------
# Plot
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(13, 7.5))

masked = np.ma.masked_invalid(
    f_matrix.values.astype(float)
)

im = ax.imshow(
    masked,
    aspect="auto",
    vmin=0,
    vmax=1
)

ax.set_xticks(np.arange(len(markers)))
ax.set_xticklabels(
    [display_names[m] for m in markers],
    rotation=35,
    ha="right"
)

ax.set_yticks(np.arange(len(genera)))
ax.set_yticklabels(genera)

ax.set_xlabel("Biomarker")
ax.set_ylabel("Fungal genus")

ax.set_title(
    "Species-level resolving performance of ITS and alternative biomarkers"
)

# --------------------------------------------------
# Cell annotations
# --------------------------------------------------

for i, genus in enumerate(genera):

    # Best alternative marker only — do not include ITS
    alternative_values = (
        f_matrix.loc[genus, markers[1:]]
        .dropna()
    )

    best_alt = (
        alternative_values.idxmax()
        if len(alternative_values)
        else None
    )

    for j, marker in enumerate(markers):

        f = f_matrix.loc[genus, marker]
        n = n_matrix.loc[genus, marker]

        if pd.isna(f):
            text = "–"
        else:
            star = "*" if marker == best_alt else ""

            text = (
                f"{f:.2f}{star}\n"
                f"n={int(n)}"
            )

        ax.text(
            j,
            i,
            text,
            ha="center",
            va="center",
            fontsize=8
        )

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Optimized genus-level F-measure")

fig.text(
    0.5,
    0.01,
    "* Highest-performing alternative biomarker within that genus. "
    "n = number of species evaluated.",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=[0, 0.04, 1, 1])

OUTPUT = OUT / "main_ITS_vs_alternative_biomarkers_heatmap.png"

fig.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print("\nF-measures:\n")
print(f_matrix.round(3).to_string())

print("\nSpecies counts:\n")
print(n_matrix.to_string())

print(f"\nSaved: {OUTPUT}")