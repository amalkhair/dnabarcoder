from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

markers = [
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
    "LNS2"
]

rows = []

for marker in markers:
    summary = Path(f"results/{marker}_focus/{marker}_indistinguishable_summary.csv")

    if summary.exists():
        df = pd.read_csv(summary)
        rows.append({
            "Marker": marker,
            "Groups": int(df["n_indistinguishable_groups"].iloc[0])
        })
    else:
        rows.append({
            "Marker": marker,
            "Groups": 0
        })

plot_df = pd.DataFrame(rows)

plt.figure(figsize=(10,5))

bars = plt.bar(
    plot_df["Marker"],
    plot_df["Groups"]
)

plt.ylabel("Number of indistinguishable species groups")
plt.xlabel("Biomarker")
plt.title("Indistinguishable species groups across fungal biomarkers")

plt.xticks(rotation=45)

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x()+bar.get_width()/2,
        height,
        f"{int(height)}",
        ha="center",
        va="bottom",
        fontsize=8
    )

plt.tight_layout()

Path("results/figures").mkdir(exist_ok=True)

plt.savefig(
    "results/figures/indistinguishable_groups_by_marker.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(plot_df)
print("\nSaved to results/figures/indistinguishable_groups_by_marker.png")