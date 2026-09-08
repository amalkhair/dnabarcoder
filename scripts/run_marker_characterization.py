from pathlib import Path
import subprocess

markers = [
    "ACT", "CAL", "LSU", "PGK", "RPB1", "RPB2",
    "SSU", "TOP1", "TUB2", "TEF1"
]

for marker in markers:
    print(f"\n=== Processing {marker} ===")

    fasta = Path(f"input/publication_scale_clean/{marker}.fasta")
    classification = Path(f"input/publication_scale_classification/{marker}.classification")
    out_dir = Path(f"results/{marker}_focus/dnabarcoder_cluster")
    clustered = out_dir / f"{marker}.clustered"

    if not fasta.exists():
        print(f"Skipping {marker}: missing FASTA")
        continue

    if not classification.exists():
        print(f"Skipping {marker}: missing classification")
        continue

    subprocess.run(
        ["python", "scripts/marker_length_stats.py", marker],
        check=True
    )

    if clustered.exists():
        print(f"Skipping clustering for {marker}: clustered file already exists")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)

        subprocess.run([
            "python", "dnabarcoder.py", "cluster",
            "-i", str(fasta),
            "-t", "1.0",
            "-rank", "species",
            "-c", str(classification),
            "-o", str(out_dir),
            "-saveindistinct", "yes"
        ], check=True)

    subprocess.run(
        ["python", "scripts/summarize_indistinct.py", marker],
        check=True
    )

print("\nAll marker characterization finished.")