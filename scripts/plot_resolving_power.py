from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path("results/final_summary/marker_threshold_summary.csv")
output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

if not input_file.exists():
    raise FileNotFoundError(f"Missing input file: {input_file}")

# Try reading with a normal header first.
df = pd.read_csv(input_file)

# Handle a CSV that may have been saved without column headers.
if "biomarker" not in df.columns:
    df = pd.read_csv(
        input_file,
        header=None,
        names=[
            "biomarker",
            "n_predicted_cutoffs",
            "mean_threshold",
            "median_threshold",
            "min_threshold",
            "max_threshold",
            "mean_f_measure",
            "median_f_measure",
            "n_sequences_used",
        ],
    )

required_columns = {
    "biomarker",
    "mean_f_measure",
    "median_f_measure",
    "n_predicted_cutoffs",
}

missing = required_columns - set(df.columns)
if missing:
    raise ValueError(
        f"Missing required columns: {sorted(missing)}\n"
        f"Available columns: {list(df.columns)}"
    )

# Convert numeric fields safely.
for column in [
    "mean_f_measure",
    "median_f_measure",
    "n_predicted_cutoffs",
]:
    df[column] = pd.to_numeric(df[column], errors="coerce")

# Remove markers for which no genus-level prediction was available.
plot_df = df[
    (df["n_predicted_cutoffs"] > 0)
    & df["mean_f_measure"].notna()
].copy()

marker_order = [
    "ITS",
    "TEF1",
    "RPB1",
    "RPB2",
    "TUB2",
    "ACT",
    "LSU",
    "SSU",
    "CAL",
    "PGK",
    "TOP1",
    "LNS2",
]

plot_df["biomarker"] = pd.Categorical(
    plot_df["biomarker"],
    categories=marker_order,
    ordered=True,
)

plot_df = plot_df.sort_values("biomarker")

fig, ax = plt.subplots(figsize=(11, 6))

bars = ax.bar(
    plot_df["biomarker"].astype(str),
    plot_df["mean_f_measure"],
)

ax.set_xlabel("Biomarker")
ax.set_ylabel("Mean F-measure")
ax.set_title("Mean resolving power across fungal biomarkers")
ax.set_ylim(0, 1.0)

ax.tick_params(axis="x", rotation=45)

for bar, value in zip(bars, plot_df["mean_f_measure"]):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.01,
        f"{value:.3f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )

fig.tight_layout()

output_file = output_dir / "mean_resolving_power.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)

# Save the underlying figure data for reproducibility.
plot_df[
    [
        "biomarker",
        "mean_f_measure",
        "median_f_measure",
        "n_predicted_cutoffs",
    ]
].to_csv(
    output_dir / "resolving_power_figure_data.csv",
    index=False,
)

print(
    plot_df[
        [
            "biomarker",
            "mean_f_measure",
            "median_f_measure",
            "n_predicted_cutoffs",
        ]
    ].to_string(index=False)
)

print(f"\nSaved figure to: {output_file}")
print(
    "Saved figure data to: "
    "results/figures/resolving_power_figure_data.csv"
)