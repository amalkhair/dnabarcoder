from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path("results/final_summary/all_marker_genus_cutoffs.csv")
output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

# Load data
df = pd.read_csv(input_file)

markers = ["RPB1", "RPB2", "TEF1"]
df = df[df["marker"].isin(markers)].copy()

# DNABarcoder "confidence" = optimized F-measure
matrix = df.pivot_table(
    index="Dataset",
    columns="marker",
    values="confidence",
    aggfunc="first"
)

comparisons = [
    ("RPB1", "RPB2"),
    ("RPB1", "TEF1"),
    ("RPB2", "TEF1")
]

summary_rows = []

rng = np.random.default_rng(42)

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False)

for ax, (m1, m2) in zip(axes, comparisons):
    sub = matrix[[m1, m2]].dropna().copy()
    sub["difference"] = sub[m1] - sub[m2]

    # Sort only for saving, not needed visually
    sub_sorted = sub.sort_values("difference").reset_index()

    # Small vertical jitter so overlapping points can be seen
    y = rng.normal(loc=0, scale=0.04, size=len(sub))

    ax.scatter(sub["difference"], y, s=35, alpha=0.8)

    # Zero line
    ax.axvline(0, linestyle="--", linewidth=1)

    # Clean up y-axis
    ax.set_yticks([])
    ax.set_ylabel("")
    ax.set_title(f"{m1} vs {m2}")

    ax.set_xlabel(f"Optimized F-measure difference ({m1} - {m2})")

    n_shared = len(sub)
    n_m1_better = (sub["difference"] > 0).sum()
    n_m2_better = (sub["difference"] < 0).sum()
    n_equal = (sub["difference"] == 0).sum()
    mean_diff = sub["difference"].mean()
    median_diff = sub["difference"].median()

    # Set symmetric-ish x-limits with small margin
    xmin = sub["difference"].min()
    xmax = sub["difference"].max()
    margin = 0.03
    ax.set_xlim(xmin - margin, xmax + margin)

    # Summary text box
    summary_text = (
        f"Shared genera: {n_shared}\n"
        f"{m1} better: {n_m1_better}\n"
        f"{m2} better: {n_m2_better}\n"
        f"Equal: {n_equal}\n"
        f"Mean diff: {mean_diff:.3f}\n"
        f"Median diff: {median_diff:.3f}"
    )

    ax.text(
        0.98, 0.95, summary_text,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9)
    )

    summary_rows.append({
        "comparison": f"{m1} vs {m2}",
        "shared_genera": n_shared,
        f"{m1}_better": n_m1_better,
        f"{m2}_better": n_m2_better,
        "equal": n_equal,
        "mean_difference": mean_diff,
        "median_difference": median_diff
    })

    # Save full genus-level differences
    sub_sorted.to_csv(
        output_dir / f"{m1}_vs_{m2}_shared_genera_differences_clean.csv",
        index=False
    )

fig.suptitle(
    "Cleaner head-to-head comparison of candidate fungal biomarkers",
    fontsize=14
)

fig.tight_layout(rect=[0, 0, 1, 0.97])

output_file = output_dir / "pairwise_marker_differences_clean.png"
fig.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close(fig)

summary_df = pd.DataFrame(summary_rows)
summary_file = output_dir / "pairwise_marker_difference_summary_clean.csv"
summary_df.to_csv(summary_file, index=False)

print(summary_df.round(4).to_string(index=False))
print(f"\nSaved figure: {output_file}")
print(f"Saved summary: {summary_file}")