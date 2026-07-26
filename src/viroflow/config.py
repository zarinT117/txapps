from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import yaml

from .models import CodingRegion, Marker, Sample


class ConfigError(ValueError):
    """Raised for a user-correctable configuration problem."""


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path} must contain a YAML mapping")
    data["_config_dir"] = str(path.resolve().parent)
    return data


def resolve_path(config: dict[str, Any], value: str | Path | None) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value)
    if not path.is_absolute():
        path = Path(config["_config_dir"]) / path
    return path.resolve()


def load_samples(config: dict[str, Any]) -> list[Sample]:
    workflow = config.get("workflow", {})
    sheet = resolve_path(config, workflow.get("samples"))
    if sheet is None:
        raise ConfigError("workflow.samples is required")
    if not sheet.exists():
        raise ConfigError(f"sample sheet does not exist: {sheet}")

    samples: list[Sample] = []
    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "r1"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ConfigError("sample sheet requires columns: sample,r1 (optional: r2,platform)")
        for row_number, row in enumerate(reader, start=2):
            name = (row.get("sample") or "").strip()
            if not name:
                raise ConfigError(f"missing sample name on row {row_number}")
            r1 = resolve_path(config, (row.get("r1") or "").strip())
            r2 = resolve_path(config, (row.get("r2") or "").strip())
            platform = (row.get("platform") or workflow.get("platform") or "illumina").lower()
            if platform not in {"illumina", "nanopore"}:
                raise ConfigError(f"unsupported platform {platform!r} for {name}")
            if r1 is None:
                raise ConfigError(f"missing r1 on row {row_number}")
            samples.append(Sample(name=name, r1=r1, r2=r2, platform=platform))
    if not samples:
        raise ConfigError("sample sheet contains no samples")
    if len({sample.name for sample in samples}) != len(samples):
        raise ConfigError("sample names must be unique")
    return samples


def parse_regions(profile: dict[str, Any]) -> list[CodingRegion]:
    regions: list[CodingRegion] = []
    for item in profile.get("coding_regions", []):
        region = CodingRegion(
            name=str(item["name"]),
            segment=str(item["segment"]),
            start=int(item["start"]),
            end=int(item["end"]),
            strand=int(item.get("strand", 1)),
        )
        if region.start < 1 or region.end < region.start or region.strand not in {-1, 1}:
            raise ConfigError(f"invalid coding region: {item}")
        regions.append(region)
    return regions


def parse_markers(profile: dict[str, Any]) -> dict[str, list[Marker]]:
    parsed: dict[str, list[Marker]] = {}
    for protein, items in profile.get("escape_markers", {}).items():
        parsed[str(protein)] = [
            Marker(
                mutation=str(item["mutation"]),
                weight=float(item.get("weight", 1.0)),
                evidence=str(item.get("evidence", "")).strip(),
            )
            for item in items
        ]
        for marker in parsed[str(protein)]:
            if not marker.evidence:
                raise ConfigError(
                    f"escape marker {protein}:{marker.mutation} requires an evidence citation"
                )
            if marker.weight < 0:
                raise ConfigError(f"escape marker weights must be non-negative: {protein}")
    return parsed

