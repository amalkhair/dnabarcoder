from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path(
    "results/final_summary/all_marker_genus_cutoffs.csv"
)

output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

markers = [
    "RPB1",
    "RPB2",
    "TEF1",
    "SSU",
    "LSU",
    "TUB2"
]

df = df[df["marker"].isin(markers)].copy()

# DNABarcoder stores the optimized F-measure
# in the column called "confidence".
data = [
    df.loc[df["marker"] == marker, "confidence"].dropna()
    for marker in markers
]

fig, ax = plt.subplots(figsize=(9, 6))

ax.boxplot(
    data,
    labels=markers,
    showfliers=True
)

ax.set_ylabel("Optimized genus-level F-measure")
ax.set_xlabel("Biomarker")
ax.set_title(
    "Distribution of species-level resolving performance\n"
    "across fungal genera"
)

ax.set_ylim(0, 1.02)

fig.tight_layout()

output_file = output_dir / "genus_fmeasure_distributions.png"
fig.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close(fig)

# Also save summary statistics
summary = (
    df.groupby("marker")["confidence"]
    .agg(
        n="count",
        mean="mean",
        median="median",
        minimum="min",
        maximum="max"
    )
)

summary["q1"] = df.groupby("marker")["confidence"].quantile(0.25)
summary["q3"] = df.groupby("marker")["confidence"].quantile(0.75)

summary = summary.loc[markers]

summary.to_csv(
    output_dir / "genus_fmeasure_distribution_summary.csv"
)

print(summary.round(4))
print(f"\nSaved: {output_file}")