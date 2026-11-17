from pathlib import Path
import pandas as pd

RESULTS = Path("results")
OUTPUT_DIR = RESULTS / "krona"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

markers = [
    "ACT",
    "CAL",
    "LSU",
    "PGK",
    "RPB1",
    "RPB2",
    "SSU",
    "TEF1",
    "TOP1",
    "TUB2",
]

rows = []

for marker in markers:

    input_file = (
        RESULTS
        / f"{marker}_focus"
        / f"{marker}_indistinguishable_groups_by_genus.csv"
    )

    if not input_file.exists():
        print(f"Missing: {input_file}")
        continue

    df = pd.read_csv(input_file)

    for _, row in df.iterrows():

        genus = str(row["genus"]).strip()
        count = int(row["n_indistinguishable_groups"])

        if not genus or count <= 0:
            continue

        rows.append({
            "count": count,
            "marker": marker,
            "genus": genus
        })


out = pd.DataFrame(rows)

# Krona ImportText format:
# count <tab> level1 <tab> level2
output_file = OUTPUT_DIR / "indistinguishable_by_marker_genus.txt"

with open(output_file, "w") as f:
    for _, row in out.iterrows():
        f.write(
            f"{row['count']}\t"
            f"{row['marker']}\t"
            f"{row['genus']}\n"
        )

out.to_csv(
    OUTPUT_DIR / "indistinguishable_by_marker_genus.csv",
    index=False
)

print(out.sort_values(
    ["marker", "count"],
    ascending=[True, False]
).to_string(index=False))

print(f"\nSaved Krona input: {output_file}")