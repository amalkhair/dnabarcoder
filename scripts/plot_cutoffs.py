import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

marker = sys.argv[1]

file = f"output/publication_scale/no_complexes/{marker}/{marker}_no_complexes.cutoffs.json.txt"

df = pd.read_csv(file, sep="\t")

outdir = Path(f"output/publication_scale/no_complexes/{marker}/figures")
outdir.mkdir(parents=True, exist_ok=True)

# Histogram
plt.figure(figsize=(8,5))
plt.hist(df["cut-off"], bins=25)
plt.xlabel("Predicted species similarity cut-off")
plt.ylabel("Number of genera")
plt.title(f"{marker} cutoff distribution")
plt.tight_layout()
plt.savefig(outdir/f"{marker}_histogram.png", dpi=300)
plt.close()

# Top genera
top = df.sort_values("sequence number", ascending=False).head(25)
top = top.sort_values("cut-off")

plt.figure(figsize=(10,7))
plt.barh(top["Dataset"], top["cut-off"])
plt.xlabel("Predicted cut-off")
plt.ylabel("Genus")
plt.title(f"{marker}: top 25 genera")
plt.tight_layout()
plt.savefig(outdir/f"{marker}_top25.png", dpi=300)
plt.close()

# Scatter
plt.figure(figsize=(8,5))
plt.scatter(df["sequence number"], df["cut-off"])
plt.xlabel("Sequence number")
plt.ylabel("Predicted cut-off")
plt.title(f"{marker}: cutoff vs dataset size")
plt.tight_layout()
plt.savefig(outdir/f"{marker}_scatter.png", dpi=300)
plt.close()

print("Done")