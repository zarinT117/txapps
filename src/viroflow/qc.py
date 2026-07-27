from __future__ import annotations

import json
import statistics
from pathlib import Path
from typing import Any

from .fasta import read_fasta


def assembly_metrics(path: Path) -> dict[str, Any]:
    records = read_fasta(path)
    lengths = sorted((len(sequence) for sequence in records.values()), reverse=True)
    total = sum(lengths)
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return {
        "contigs": len(lengths),
        "total_bases": total,
        "largest_contig": lengths[0],
        "n50": n50,
        "gc_fraction": round(
            sum(sequence.count("G") + sequence.count("C") for sequence in records.values())
            / total,
            6,
        )
        if total
        else 0.0,
    }


def coverage_metrics(path: Path, thresholds: tuple[int, ...] = (1, 10, 100)) -> dict[str, Any]:
    depths: list[int] = []
    contigs: dict[str, list[int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                raise ValueError(f"invalid depth row {line_number} in {path}")
            contig, depth = parts[0], int(parts[2])
            depths.append(depth)
            contigs.setdefault(contig, []).append(depth)
    if not depths:
        raise ValueError(f"depth file contains no positions: {path}")
    mean_depth = statistics.fmean(depths)
    result: dict[str, Any] = {
        "positions": len(depths),
        "mean_depth": round(mean_depth, 3),
        "median_depth": round(statistics.median(depths), 3),
        "min_depth": min(depths),
        "max_depth": max(depths),
        "uniformity_fraction_0_2x_to_2x_mean": round(
            sum(0.2 * mean_depth <= depth <= 2 * mean_depth for depth in depths) / len(depths),
            6,
        )
        if mean_depth
        else 0.0,
        "contigs": {},
    }
    for threshold in thresholds:
        result[f"breadth_at_{threshold}x"] = round(
            sum(depth >= threshold for depth in depths) / len(depths), 6
        )
    for contig, values in contigs.items():
        result["contigs"][contig] = {
            "positions": len(values),
            "mean_depth": round(statistics.fmean(values), 3),
            "breadth_at_10x": round(sum(value >= 10 for value in values) / len(values), 6),
        }
    return result


def write_qc_summary(
    output_path: Path,
    assembly_path: Path,
    depth_path: Path,
    additional: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = {
        "assembly": assembly_metrics(assembly_path),
        "coverage": coverage_metrics(depth_path),
        **(additional or {}),
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload

