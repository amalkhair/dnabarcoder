from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


BASE = Path("results/barcode_gap")

DATASETS = {
    "ITS": {
        "sim": BASE / "Fusarium_ITS/Fusarium_ITS_matched_balanced.sim",
        "fasta": Path(
            "results/phylo_validation/fusarium_ITS_vs_CaM/"
            "Fusarium_ITS_matched_balanced.fasta"
        ),
    },
    "CaM": {
        "sim": BASE / "Fusarium_CaM/Fusarium_CaM_matched_balanced.sim",
        "fasta": Path(
            "results/phylo_validation/fusarium_ITS_vs_CaM/"
            "Fusarium_CaM_matched_balanced.fasta"
        ),
    },
}

OUT = BASE / "Fusarium_ITS_vs_CaM"
OUT.mkdir(parents=True, exist_ok=True)


def species_from_id(seq_id):
    """
    Example:
    Fusarium_acuminatum|ITS|rep1|PP959301.1
    -> Fusarium_acuminatum
    """
    return seq_id.split("|")[0]


def read_fasta_ids(path):
    ids = []

    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])

    return ids


def read_similarity_file(path):
    rows = []

    with open(path) as handle:
        for line in handle:
            parts = line.strip().split()

            if len(parts) != 3:
                continue

            seq1, seq2, similarity = parts

            rows.append({
                "seq1": seq1,
                "seq2": seq2,
                "similarity": float(similarity),
            })

    df = pd.DataFrame(rows)

    # Remove self comparisons
    df = df[df["seq1"] != df["seq2"]].copy()

    df["species1"] = df["seq1"].map(species_from_id)
    df["species2"] = df["seq2"].map(species_from_id)

    return df


def analyse_marker(marker, sim_path, fasta_path):

    print(f"\n{'=' * 70}")
    print(marker)
    print("=" * 70)

    sim = read_similarity_file(sim_path)

    fasta_ids = read_fasta_ids(fasta_path)

    fasta_species = {}

    for seq_id in fasta_ids:
        species = species_from_id(seq_id)
        fasta_species.setdefault(species, []).append(seq_id)

    results = []

    for species in sorted(fasta_species):

        # ------------------------------------------
        # INTRASPECIFIC COMPARISONS
        # ------------------------------------------

        intra = sim[
            (sim["species1"] == species)
            & (sim["species2"] == species)
        ].copy()

        # Number expected from n sequences:
        # n choose 2
        n_seq = len(fasta_species[species])
        expected_intra_pairs = n_seq * (n_seq - 1) // 2

        # Remove possible duplicated A-B / B-A pairs
        intra_pairs = {}

        for _, row in intra.iterrows():

            pair = tuple(sorted([
                row["seq1"],
                row["seq2"],
            ]))

            intra_pairs.setdefault(pair, []).append(
                row["similarity"]
            )

        intra_values = [
            sum(values) / len(values)
            for values in intra_pairs.values()
        ]

        if not intra_values:
            print(
                f"WARNING: no within-species comparisons for {species}"
            )
            continue

        # Lowest within-species similarity corresponds
        # to maximum within-species distance.
        min_intra_similarity = min(intra_values)

        max_intra_distance = 1 - min_intra_similarity


        # ------------------------------------------
        # INTERSPECIFIC COMPARISONS
        # ------------------------------------------

        inter = sim[
            (
                (sim["species1"] == species)
                & (sim["species2"] != species)
            )
            |
            (
                (sim["species2"] == species)
                & (sim["species1"] != species)
            )
        ].copy()

        if inter.empty:
            print(
                f"WARNING: no between-species comparisons for {species}"
            )
            continue

        # Find the most similar sequence belonging
        # to another species.
        nearest_row = inter.loc[
            inter["similarity"].idxmax()
        ]

        max_inter_similarity = nearest_row["similarity"]

        if nearest_row["species1"] == species:
            nearest_species = nearest_row["species2"]
        else:
            nearest_species = nearest_row["species1"]

        nearest_inter_distance = 1 - max_inter_similarity


        # ------------------------------------------
        # BARCODE GAP
        # ------------------------------------------

        barcode_gap = (
            nearest_inter_distance
            - max_intra_distance
        )

        # Equivalent:
        # min_intra_similarity - max_inter_similarity

        results.append({
            "marker": marker,
            "species": species,
            "n_sequences": n_seq,

            "expected_intra_pairs": expected_intra_pairs,
            "observed_intra_pairs": len(intra_pairs),

            "min_intraspecific_similarity":
                min_intra_similarity,

            "max_intraspecific_distance":
                max_intra_distance,

            "nearest_other_species":
                nearest_species,

            "nearest_interspecific_similarity":
                max_inter_similarity,

            "nearest_interspecific_distance":
                nearest_inter_distance,

            "barcode_gap":
                barcode_gap,

            "positive_barcode_gap":
                barcode_gap > 0,
        })

    result = pd.DataFrame(results)

    result.to_csv(
        OUT / f"Fusarium_{marker}_barcode_gap_by_species.csv",
        index=False
    )

    return result


# ==================================================
# RUN BOTH MARKERS
# ==================================================

all_results = []

for marker, paths in DATASETS.items():

    result = analyse_marker(
        marker,
        paths["sim"],
        paths["fasta"],
    )

    all_results.append(result)


combined = pd.concat(
    all_results,
    ignore_index=True
)

combined.to_csv(
    OUT / "Fusarium_ITS_CaM_barcode_gap_all.csv",
    index=False
)


# ==================================================
# SUMMARY
# ==================================================

summary_rows = []

for marker, group in combined.groupby("marker"):

    positive = group["positive_barcode_gap"].sum()
    total = len(group)

    summary_rows.append({
        "marker": marker,
        "species_evaluated": total,
        "positive_gap_species": positive,
        "positive_gap_percent":
            positive / total * 100,
        "median_barcode_gap":
            group["barcode_gap"].median(),
        "mean_barcode_gap":
            group["barcode_gap"].mean(),
        "median_max_intra_distance":
            group["max_intraspecific_distance"].median(),
        "median_nearest_inter_distance":
            group["nearest_interspecific_distance"].median(),
    })


summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT / "Fusarium_ITS_CaM_barcode_gap_summary.csv",
    index=False
)

print("\nBARCODE GAP SUMMARY")
print("=" * 70)
print(summary.to_string(index=False))


# ==================================================
# MATCH ITS AND CaM SPECIES
# ==================================================

its = combined[
    combined["marker"] == "ITS"
][
    [
        "species",
        "barcode_gap",
        "positive_barcode_gap",
    ]
].rename(
    columns={
        "barcode_gap": "ITS_gap",
        "positive_barcode_gap": "ITS_positive",
    }
)

cam = combined[
    combined["marker"] == "CaM"
][
    [
        "species",
        "barcode_gap",
        "positive_barcode_gap",
    ]
].rename(
    columns={
        "barcode_gap": "CaM_gap",
        "positive_barcode_gap": "CaM_positive",
    }
)

paired = its.merge(
    cam,
    on="species",
    how="inner"
)

paired["gap_difference_CaM_minus_ITS"] = (
    paired["CaM_gap"] - paired["ITS_gap"]
)

paired["CaM_rescued"] = (
    (~paired["ITS_positive"])
    & paired["CaM_positive"]
)

paired["ITS_only"] = (
    paired["ITS_positive"]
    & (~paired["CaM_positive"])
)

paired.to_csv(
    OUT / "Fusarium_ITS_vs_CaM_barcode_gap_paired.csv",
    index=False
)


print("\nPAIRED COMPARISON")
print("=" * 70)

print(f"Shared species: {len(paired)}")

print(
    "CaM larger barcode gap:",
    (paired["gap_difference_CaM_minus_ITS"] > 0).sum()
)

print(
    "ITS larger barcode gap:",
    (paired["gap_difference_CaM_minus_ITS"] < 0).sum()
)

print(
    "CaM rescued ITS-negative species:",
    paired["CaM_rescued"].sum()
)

print(
    "ITS positive but CaM negative:",
    paired["ITS_only"].sum()
)


# ==================================================
# FIGURE 1 — DISTRIBUTION OF SPECIES BARCODE GAPS
# ==================================================

fig, ax = plt.subplots(figsize=(7, 6))

its_values = combined.loc[
    combined["marker"] == "ITS",
    "barcode_gap"
]

cam_values = combined.loc[
    combined["marker"] == "CaM",
    "barcode_gap"
]

ax.boxplot(
    [its_values, cam_values],
    tick_labels=["ITS", "CaM"],
    showmeans=True
)

ax.axhline(
    0,
    linestyle="--",
    linewidth=1
)

ax.set_ylabel("Species-level barcode gap")
ax.set_title(
    "Barcode-gap comparison for matched Fusarium species"
)

ax.text(
    0.02,
    0.02,
    "Positive values indicate separation between\n"
    "intra- and interspecific variation.",
    transform=ax.transAxes,
    fontsize=9,
    va="bottom"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    OUT / "Fusarium_ITS_vs_CaM_barcode_gap_boxplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ==================================================
# FIGURE 2 — PERCENTAGE WITH POSITIVE GAP
# ==================================================

plot_summary = summary.sort_values("marker")

fig, ax = plt.subplots(figsize=(6, 5))

bars = ax.bar(
    plot_summary["marker"],
    plot_summary["positive_gap_percent"]
)

for bar, value in zip(
    bars,
    plot_summary["positive_gap_percent"]
):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + 1,
        f"{value:.1f}%",
        ha="center",
        fontsize=11
    )

ax.set_ylim(0, 100)

ax.set_ylabel(
    "Species with a positive barcode gap (%)"
)

ax.set_title(
    "Species-level barcode separation in Fusarium"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    OUT / "Fusarium_ITS_vs_CaM_positive_gap_percent.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print("\nSaved results to:")
print(OUT)