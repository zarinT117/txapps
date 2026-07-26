from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


def read_fasta(path: Path) -> dict[str, str]:
    records = {
        record.id: str(record.seq).upper().replace("U", "T")
        for record in SeqIO.parse(str(path), "fasta")
    }
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records


def write_fasta(records: dict[str, str], path: Path, width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), width):
                handle.write(sequence[offset : offset + width] + "\n")

