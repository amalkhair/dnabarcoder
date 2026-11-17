from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path("results/final_summary/marker_threshold_summary.csv")
output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

# Include ITS as the reference marker.
# Exclude markers with very limited coverage from this figure;
# they remain reported in the full summary table.
df = df[
    (~df["marker"].isin(["LNS2", "ACT", "CAL", "PGK", "TOP1"]))
    & (df["n_cutoffs"] > 0)
].copy()

# DNABarcoder stores the optimized F-measure
# in the output column called "confidence".
df["mean_f_measure"] = df["mean_confidence"]

fig, ax = plt.subplots(figsize=(9, 6))

ax.scatter(
    df["n_cutoffs"],
    df["mean_f_measure"],
    s=90
)

for _, row in df.iterrows():
    ax.annotate(
        row["marker"],
        (row["n_cutoffs"], row["mean_f_measure"]),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=9
    )

ax.set_xlabel("Number of genera with predicted species-level cutoff")
ax.set_ylabel("Mean optimized F-measure")

ax.set_title(
    "Taxonomic coverage vs species-level resolving performance"
)

ax.set_ylim(0.70, 1.0)

fig.tight_layout()

output_file = output_dir / "coverage_vs_resolving_performancev2.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

# Save the exact data used in the figure
df[
    ["marker", "n_cutoffs", "mean_f_measure"]
].sort_values(
    "n_cutoffs",
    ascending=False
).to_csv(
    output_dir / "coverage_vs_resolving_performance_data.csv",
    index=False
)

print(
    df[
        ["marker", "n_cutoffs", "mean_f_measure"]
    ].sort_values(
        "n_cutoffs",
        ascending=False
    ).to_string(index=False)
)

print(f"\nSaved: {output_file}")