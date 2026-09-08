import json
from pathlib import Path
from collections import Counter
import pandas as pd

input_file = Path("results/TEF1_focus/dnabarcoder_cluster/TEF1.indistinct.species")
out_dir = Path("results/TEF1_focus")
out_dir.mkdir(parents=True, exist_ok=True)

with open(input_file, "r", encoding="utf-8") as f:
    data = json.load(f)

# Remove duplicate groups
unique_groups = {}
for species, group in data.items():
    group_sorted = tuple(sorted(group))
    unique_groups[group_sorted] = group_sorted

rows = []
for group in unique_groups.values():
    genera = [name.split()[0] for name in group if name.strip()]
    most_common_genus = Counter(genera).most_common(1)[0][0] if genera else "Unknown"

    rows.append({
        "representative_genus": most_common_genus,
        "n_species_in_group": len(group),
        "species_group": "; ".join(group)
    })

df = pd.DataFrame(rows)
df = df.sort_values("n_species_in_group", ascending=False)

summary = pd.DataFrame({
    "biomarker": ["TEF1"],
    "n_indistinguishable_groups": [len(df)],
    "mean_species_per_group": [round(df["n_species_in_group"].mean(), 2)],
    "max_species_in_group": [df["n_species_in_group"].max()]
})

genus_counts = df["representative_genus"].value_counts().reset_index()
genus_counts.columns = ["genus", "n_indistinguishable_groups"]

df.to_csv(out_dir / "TEF1_indistinguishable_species_groups.csv", index=False)
summary.to_csv(out_dir / "TEF1_indistinguishable_summary.csv", index=False)
genus_counts.to_csv(out_dir / "TEF1_indistinguishable_groups_by_genus.csv", index=False)

print(summary)
print("\nTop genera:")
print(genus_counts.head(15))
print(f"\nSaved outputs to: {out_dir}")