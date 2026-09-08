from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

markers = [
    "ITS","TEF1","RPB1","RPB2","TUB2",
    "ACT","LSU","SSU","CAL","PGK","TOP1","LNS2"
]

dfs = []

for marker in markers:
    file = Path(f"results/{marker}_focus/{marker}_indistinguishable_groups_by_genus.csv")

    if not file.exists():
        continue

    df = pd.read_csv(file)
    df["marker"] = marker
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

# make sure counts are numeric
all_df["n_indistinguishable_groups"] = pd.to_numeric(
    all_df["n_indistinguishable_groups"],
    errors="coerce"
).fillna(0)

# top 20 genera overall
top_genera = (
    all_df.groupby("genus")["n_indistinguishable_groups"]
    .sum()
    .sort_values(ascending=False)
    .head(20)
    .index
)

heatmap = (
    all_df[all_df["genus"].isin(top_genera)]
    .pivot_table(
        index="genus",
        columns="marker",
        values="n_indistinguishable_groups",
        aggfunc="sum",
        fill_value=0
    )
)

# keep biomarker order consistent
heatmap = heatmap.reindex(columns=[m for m in markers if m in heatmap.columns])

# sort genera by total counts
heatmap["TOTAL"] = heatmap.sum(axis=1)
heatmap = heatmap.sort_values("TOTAL", ascending=False)
heatmap = heatmap.drop(columns="TOTAL")

# convert to numeric numpy array
heatmap = heatmap.astype(float)

plt.figure(figsize=(11,8))

im = plt.imshow(
    heatmap.values,
    cmap="viridis",
    aspect="auto"
)

plt.colorbar(im, label="Indistinguishable groups")

plt.xticks(
    range(len(heatmap.columns)),
    heatmap.columns,
    rotation=45,
    ha="right"
)

plt.yticks(
    range(len(heatmap.index)),
    heatmap.index,
    fontstyle="italic"
)

plt.xlabel("Biomarker")
plt.ylabel("Genus")
plt.title("Top genera associated with indistinguishable species groups")

plt.tight_layout()

Path("results/figures").mkdir(exist_ok=True)

plt.savefig(
    "results/figures/genus_indistinguishable_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:")
print("results/figures/genus_indistinguishable_heatmap.png")