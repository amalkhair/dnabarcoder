from pathlib import Path
from Bio import SeqIO

in_dir = Path("input/publication_scale")
fasta_out_dir = Path("input/publication_scale_clean")
class_out_dir = Path("input/publication_scale_classification")

fasta_out_dir.mkdir(parents=True, exist_ok=True)
class_out_dir.mkdir(parents=True, exist_ok=True)

for fasta in sorted(in_dir.glob("*_no_complexes.fasta")):
    marker = fasta.name.replace("_no_complexes.fasta", "")

    out_fasta = fasta_out_dir / f"{marker}.fasta"
    out_class = class_out_dir / f"{marker}.classification"

    records_out = []

    with open(out_class, "w") as c:
        c.write("id\tgenus\tspecies\n")

        for record in SeqIO.parse(fasta, "fasta"):
            seq_id = record.id
            parts = record.description.split()

            if len(parts) < 3:
                continue

            genus = parts[1]
            species_epithet = parts[2]
            species = f"{genus} {species_epithet}"

            record.id = seq_id
            record.name = seq_id
            record.description = ""

            records_out.append(record)
            c.write(f"{seq_id}\t{genus}\t{species}\n")

    SeqIO.write(records_out, out_fasta, "fasta")

    print(marker, len(records_out))
