import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =========================================================
# OUTPUT
# =========================================================

OUT_DIR = Path("results/visualisation_analysis/section_3_1")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =========================================================
# DATA
# Same subset is used for all three visualisations
# =========================================================

df = pd.DataFrame(
    {
        "ITS":  [0.724, 0.668, 0.567, 0.604],
        "RPB1": [0.961, 0.925, 0.732, np.nan],
        "RPB2": [0.840, 0.833, 0.682, 0.772],
        "TEF1": [0.890, 0.913, 0.722, 0.669],
        "TUB2": [0.832, 0.778, 0.625, 0.782],
        "CaM":  [0.857, 0.810, 0.740, 0.719],
    },
    index=[
        "Aspergillus",
        "Penicillium",
        "Fusarium",
        "Trichoderma",
    ],
)

df.index.name = "Genus"

# Save exact dataset used
df.to_csv(OUT_DIR / "figure_3_1_data.csv")


# =========================================================
# FIGURE 3.1A — TABLE
# =========================================================

fig, ax = plt.subplots(figsize=(9, 3.2))
ax.axis("off")

cell_text = [
    [
        f"{value:.3f}" if pd.notna(value) else "N/A"
        for value in row
    ]
    for row in df.to_numpy()
]

table = ax.table(
    cellText=cell_text,
    rowLabels=df.index,
    colLabels=df.columns,
    cellLoc="center",
    loc="center",
)

table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 1.6)

ax.set_title(
    "Genus-level biomarker F-measures — Table",
    pad=15,
)

plt.tight_layout()

plt.savefig(
    OUT_DIR / "figure_3_1A_table.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# =========================================================
# FIGURE 3.1B — GROUPED BAR CHART
# =========================================================

fig, ax = plt.subplots(figsize=(11, 6))

x = np.arange(len(df.index))

n_markers = len(df.columns)
bar_width = 0.12

for i, marker in enumerate(df.columns):

    offset = (
        i - (n_markers - 1) / 2
    ) * bar_width

    values = df[marker].to_numpy()

    ax.bar(
        x + offset,
        values,
        width=bar_width,
        label=marker,
    )

ax.set_xticks(x)
ax.set_xticklabels(df.index)

ax.set_xlabel("Genus")
ax.set_ylabel("F-measure")

ax.set_ylim(0, 1)

ax.set_title(
    "Genus-level biomarker F-measures — Grouped Bar Chart"
)

ax.legend(
    title="Biomarker",
    ncol=3,
    frameon=False,
)

plt.tight_layout()

plt.savefig(
    OUT_DIR / "figure_3_1B_grouped_bar.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# =========================================================
# FIGURE 3.1C — ANNOTATED HEATMAP
# =========================================================

fig, ax = plt.subplots(figsize=(9, 5))

values = df.to_numpy()

# Mask missing values
masked_values = np.ma.masked_invalid(values)

heatmap = ax.imshow(
    masked_values,
    aspect="auto",
    vmin=0,
    vmax=1,
)

# Axis labels
ax.set_xticks(np.arange(len(df.columns)))
ax.set_xticklabels(df.columns)

ax.set_yticks(np.arange(len(df.index)))
ax.set_yticklabels(df.index)

ax.set_xlabel("Biomarker")
ax.set_ylabel("Genus")

ax.set_title(
    "Genus-level biomarker F-measures — Heatmap"
)

# Add values to cells
for row in range(df.shape[0]):
    for col in range(df.shape[1]):

        value = df.iloc[row, col]

        if pd.isna(value):
            label = "N/A"
        else:
            label = f"{value:.3f}"

        ax.text(
            col,
            row,
            label,
            ha="center",
            va="center",
        )

# Color bar
colorbar = fig.colorbar(
    heatmap,
    ax=ax,
)

colorbar.set_label("F-measure")

plt.tight_layout()

plt.savefig(
    OUT_DIR / "figure_3_1C_heatmap.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()


print("Done.")
print(f"Figures saved to: {OUT_DIR}")