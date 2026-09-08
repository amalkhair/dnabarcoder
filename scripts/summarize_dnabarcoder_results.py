from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

markers = ["ACT","CAL","ITS","LNS2","LSU","PGK","RPB1","RPB2","SSU","TEF1","TOP1","TUB2"]

base = Path("output/publication_scale/no_complexes")
summary_dir = Path("results/final_summary")
fig_dir = Path("results/figures")

summary_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

all_rows = []
summary_rows = []

for marker in markers:
    file = base / marker / f"{marker}_no_complexes.cutoffs.json.txt"

    if not file.exists():
        summary_rows.append({
            "marker": marker,
            "n_cutoffs": 0,
            "mean_cutoff": None,
            "median_cutoff": None,
            "min_cutoff": None,
            "max_cutoff": None,
            "mean_confidence": None,
            "median_confidence": None,
            "total_sequences_in_thresholded_genera": 0,
        })
        continue

    df = pd.read_csv(file, sep="\t")

    if df.empty:
        n = 0
        summary_rows.append({
            "marker": marker,
            "n_cutoffs": 0,
            "mean_cutoff": None,
            "median_cutoff": None,
            "min_cutoff": None,
            "max_cutoff": None,
            "mean_confidence": None,
            "median_confidence": None,
            "total_sequences_in_thresholded_genera": 0,
        })
        continue

    df["marker"] = marker
    all_rows.append(df)

    summary_rows.append({
        "marker": marker,
        "n_cutoffs": len(df),
        "mean_cutoff": df["cut-off"].mean(),
        "median_cutoff": df["cut-off"].median(),
        "min_cutoff": df["cut-off"].min(),
        "max_cutoff": df["cut-off"].max(),
        "mean_confidence": df["confidence"].mean(),
        "median_confidence": df["confidence"].median(),
        "total_sequences_in_thresholded_genera": df["sequence number"].sum(),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(summary_dir / "marker_threshold_summary.csv", index=False)

if all_rows:
    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(summary_dir / "all_marker_genus_cutoffs.csv", index=False)

    # Figure 1: number of cutoffs
    plot_df = summary.sort_values("n_cutoffs", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.bar(plot_df["marker"], plot_df["n_cutoffs"])
    plt.xlabel("Marker")
    plt.ylabel("Number of predicted genus-level cutoffs")
    plt.title("DNABarcoder predicted cutoffs per marker")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(fig_dir / "marker_cutoff_counts.png", dpi=300)
    plt.close()

    # Figure 2: threshold distributions
    dist_df = combined[combined["marker"].isin(
        ["ITS","LSU","SSU","TEF1","RPB1","RPB2","TUB2","ACT","CAL"]
    )]

    order = ["ITS","LSU","SSU","TEF1","RPB1","RPB2","TUB2","ACT","CAL"]
    data = [dist_df.loc[dist_df["marker"] == m, "cut-off"].dropna() for m in order]

    plt.figure(figsize=(11, 6))
    plt.boxplot(data, labels=order, showfliers=False)
    plt.xlabel("Marker")
    plt.ylabel("Predicted species similarity cutoff")
    plt.title("Distribution of predicted cutoffs across fungal biomarkers")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(fig_dir / "marker_cutoff_distributions_boxplot.png", dpi=300)
    plt.close()

print("Saved:")
print(summary_dir / "marker_threshold_summary.csv")
print(summary_dir / "all_marker_genus_cutoffs.csv")
print(fig_dir / "marker_cutoff_counts.png")
print(fig_dir / "marker_cutoff_distributions_boxplot.png")