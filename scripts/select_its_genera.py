from pathlib import Path
import pandas as pd
import numpy as np

INPUT = Path(
    "results/final_summary/all_marker_genus_cutoffs.csv"
)

OUTPUT_DIR = Path("results/genus_selection")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

# --------------------------------------------------
# Keep ITS only
# --------------------------------------------------

its = df[
    df["marker"].astype(str).str.upper() == "ITS"
].copy()

# Make sure relevant columns are numeric
its["cut-off"] = pd.to_numeric(
    its["cut-off"], errors="coerce"
)

its["sequence number"] = pd.to_numeric(
    its["sequence number"], errors="coerce"
)

its["group number"] = pd.to_numeric(
    its["group number"], errors="coerce"
)

its["confidence"] = pd.to_numeric(
    its["confidence"], errors="coerce"
)

# --------------------------------------------------
# Duong's selection rule
#
# 1. ITS optimized threshold = 1.000
# 2. At least 1000 sequences
# --------------------------------------------------

threshold_100 = its[
    np.isclose(its["cut-off"], 1.000)
].copy()

selected = threshold_100[
    threshold_100["sequence number"] >= 1000
].copy()

excluded_low_n = threshold_100[
    threshold_100["sequence number"] < 1000
].copy()

# Sort most strongly represented genera first
selected = selected.sort_values(
    "sequence number",
    ascending=False
)

excluded_low_n = excluded_low_n.sort_values(
    "sequence number",
    ascending=False
)

# --------------------------------------------------
# Save
# --------------------------------------------------

columns = [
    "Dataset",
    "cut-off",
    "sequence number",
    "group number",
    "confidence"
]

selected[columns].to_csv(
    OUTPUT_DIR / "ITS_selected_genera.csv",
    index=False
)

excluded_low_n[columns].to_csv(
    OUTPUT_DIR / "ITS_threshold_1_excluded_under_1000.csv",
    index=False
)

threshold_100[columns].to_csv(
    OUTPUT_DIR / "ITS_all_threshold_1_genera.csv",
    index=False
)

# --------------------------------------------------
# Print results
# --------------------------------------------------

print("\n" + "=" * 75)
print("ITS GENUS SELECTION")
print("=" * 75)

print(f"\nTotal ITS genera with predicted cutoff: {len(its)}")
print(
    "ITS genera with cutoff = 1.000:       "
    f"{len(threshold_100)}"
)
print(
    "Cutoff = 1.000 AND >=1000 sequences:  "
    f"{len(selected)}"
)

print("\nSELECTED GENERA")
print("-" * 75)

if selected.empty:
    print("No genera passed both criteria.")
else:
    print(
        selected[columns]
        .to_string(index=False)
    )

print("\nEXCLUDED BECAUSE <1000 SEQUENCES")
print("-" * 75)

if excluded_low_n.empty:
    print("None")
else:
    print(
        excluded_low_n[columns]
        .to_string(index=False)
    )

print("\nSaved:")
print(
    OUTPUT_DIR / "ITS_selected_genera.csv"
)
print(
    OUTPUT_DIR /
    "ITS_threshold_1_excluded_under_1000.csv"
)