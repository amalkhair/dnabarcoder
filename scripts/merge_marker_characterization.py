from pathlib import Path
import pandas as pd

results_dir = Path("results")
rows = []

for summary_file in sorted(results_dir.glob("*_focus/*_sequence_length_summary.csv")):
    marker = summary_file.name.replace("_sequence_length_summary.csv", "")
    indistinct_file = results_dir / f"{marker}_focus" / f"{marker}_indistinguishable_summary.csv"

    length_df = pd.read_csv(summary_file)

    if indistinct_file.exists():
        indistinct_df = pd.read_csv(indistinct_file)
    else:
        indistinct_df = pd.DataFrame([{
            "biomarker": marker,
            "n_indistinguishable_groups": 0,
            "mean_species_per_group": 0,
            "max_species_in_group": 0
        }])

    merged = pd.merge(length_df, indistinct_df, on="biomarker", how="left")
    rows.append(merged)

master = pd.concat(rows, ignore_index=True)
master = master.sort_values("biomarker")

out_dir = Path("results/final_summary")
out_dir.mkdir(parents=True, exist_ok=True)

master.to_csv(out_dir / "marker_characterization_summary.csv", index=False)

print(master)
print(f"\nSaved to: {out_dir / 'marker_characterization_summary.csv'}")