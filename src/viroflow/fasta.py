from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id in records:
            raise ValueError(f"duplicate FASTA record ID {record.id!r} in {path}")
        sequence = str(record.seq).upper().replace("U", "T")
        invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN-."))
        if invalid:
            raise ValueError(
                f"invalid nucleotide symbol(s) in {path}, record {record.id}: {''.join(invalid)}"
            )
        records[record.id] = sequence
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records


def write_fasta(records: dict[str, str], path: Path, width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), width):
                handle.write(sequence[offset : offset + width] + "\n")
