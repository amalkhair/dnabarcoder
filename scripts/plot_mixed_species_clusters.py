import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

data = {
    "marker": [
        "ITS","TEF1","RPB1","RPB2",
        "TUB2","ACT","CAL","LSU",
        "SSU","PGK","TOP1","LNS2"
    ],
    "mixed_clusters": [
        3525,2067,488,1768,
        3049,363,1240,3284,
        3024,7,4,1
    ]
}

df = pd.DataFrame(data)

df = df.sort_values("mixed_clusters", ascending=False)

Path("results/figures").mkdir(parents=True, exist_ok=True)

plt.figure(figsize=(11,6))

plt.bar(df["marker"], df["mixed_clusters"])

plt.xlabel("Biomarker")
plt.ylabel("Number of mixed-species clusters")
plt.title("Mixed-species clusters identified at 100% sequence identity")

for i, v in enumerate(df["mixed_clusters"]):
    plt.text(i, v + 40, str(v), ha='center', fontsize=8)

plt.tight_layout()

plt.savefig(
    "results/figures/mixed_species_clusters_by_marker.png",
    dpi=300
)

plt.show()