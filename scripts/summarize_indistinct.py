import json
from pathlib import Path
from collections import Counter
import sys
import pandas as pd

if len(sys.argv) != 2:
    print("Usage: python scripts/summarize_indistinct.py MARKER")
    sys.exit(1)

marker = sys.argv[1].upper()

input_file = Path(f"results/{marker}_focus/dnabarcoder_cluster/{marker}.indistinct.species")
out_dir = Path(f"results/{marker}_focus")
out_dir.mkdir(parents=True, exist_ok=True)

summary_path = out_dir / f"{marker}_indistinguishable_summary.csv"
groups_path = out_dir / f"{marker}_indistinguishable_species_groups.csv"
genus_path = out_dir / f"{marker}_indistinguishable_groups_by_genus.csv"

if not input_file.exists():
    print(f"No indistinguishable species file found for {marker}. Assuming 0 groups.")

    summary = pd.DataFrame({
        "biomarker": [marker],
        "n_indistinguishable_groups": [0],
        "mean_species_per_group": [0],
        "max_species_in_group": [0]
    })

    empty_groups = pd.DataFrame(columns=[
        "representative_genus",
        "n_species_in_group",
        "species_group"
    ])

    empty_genus = pd.DataFrame(columns=[
        "genus",
        "n_indistinguishable_groups"
    ])

    summary.to_csv(summary_path, index=False)
    empty_groups.to_csv(groups_path, index=False)
    empty_genus.to_csv(genus_path, index=False)

    print(summary)
    print(f"\nSaved empty outputs to: {out_dir}")
    sys.exit(0)

with open(input_file, "r", encoding="utf-8") as f:
    data = json.load(f)

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

if df.empty:
    summary = pd.DataFrame({
        "biomarker": [marker],
        "n_indistinguishable_groups": [0],
        "mean_species_per_group": [0],
        "max_species_in_group": [0]
    })

    genus_counts = pd.DataFrame(columns=["genus", "n_indistinguishable_groups"])
else:
    df = df.sort_values("n_species_in_group", ascending=False)

    summary = pd.DataFrame({
        "biomarker": [marker],
        "n_indistinguishable_groups": [len(df)],
        "mean_species_per_group": [round(df["n_species_in_group"].mean(), 2)],
        "max_species_in_group": [df["n_species_in_group"].max()]
    })

    genus_counts = df["representative_genus"].value_counts().reset_index()
    genus_counts.columns = ["genus", "n_indistinguishable_groups"]

df.to_csv(groups_path, index=False)
summary.to_csv(summary_path, index=False)
genus_counts.to_csv(genus_path, index=False)

print(summary)
print("\nTop genera:")
print(genus_counts.head(15))
print(f"\nSaved outputs to: {out_dir}")