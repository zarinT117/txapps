from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import yaml

from .models import CodingRegion, Marker, Sample


class ConfigError(ValueError):
    """Raised for a user-correctable configuration problem."""


SAMPLE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


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
            if not SAMPLE_NAME.fullmatch(name):
                raise ConfigError(
                    f"invalid sample name {name!r} on row {row_number}; use letters, numbers, "
                    "periods, underscores, or hyphens (maximum 100 characters)"
                )
            r1 = resolve_path(config, (row.get("r1") or "").strip())
            r2 = resolve_path(config, (row.get("r2") or "").strip())
            platform = (row.get("platform") or workflow.get("platform") or "illumina").lower()
            if platform not in {"illumina", "nanopore"}:
                raise ConfigError(f"unsupported platform {platform!r} for {name}")
            if r1 is None:
                raise ConfigError(f"missing r1 on row {row_number}")
            if not str(r1).lower().endswith(FASTQ_SUFFIXES):
                raise ConfigError(f"r1 must be FASTQ or gzipped FASTQ for {name}: {r1}")
            if r2 and not str(r2).lower().endswith(FASTQ_SUFFIXES):
                raise ConfigError(f"r2 must be FASTQ or gzipped FASTQ for {name}: {r2}")
            if platform == "nanopore" and r2:
                raise ConfigError(f"Nanopore sample {name} must not define r2")
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
        if (region.end - region.start + 1) % 3:
            raise ConfigError(f"coding region length must be divisible by 3: {region.name}")
        regions.append(region)
    if len({region.name for region in regions}) != len(regions):
        raise ConfigError("coding region names must be unique")
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


def validate_config(config: dict[str, Any], check_files: bool = True) -> list[str]:
    """Validate workflow settings and return non-fatal advisory messages."""
    workflow = config.get("workflow")
    if not isinstance(workflow, dict):
        raise ConfigError("workflow must be a YAML mapping")
    samples = load_samples(config)
    reference = resolve_path(config, workflow.get("reference_fasta"))
    if reference is None:
        raise ConfigError("workflow.reference_fasta is required")
    optional_files = {
        "host_reference_fasta": resolve_path(config, workflow.get("host_reference_fasta")),
        "primer_bed": resolve_path(config, workflow.get("primer_bed")),
        "nextclade_dataset": resolve_path(config, workflow.get("nextclade_dataset")),
    }
    if check_files:
        required_paths = [reference, *(sample.r1 for sample in samples)]
        required_paths.extend(sample.r2 for sample in samples if sample.r2)
        missing = [str(path) for path in required_paths if not path.exists()]
        missing.extend(
            f"{name}={path}" for name, path in optional_files.items() if path and not path.exists()
        )
        if missing:
            raise ConfigError("required input path(s) do not exist: " + ", ".join(missing))

    numeric_rules = {
        "threads": (1, 1024, int),
        "min_read_length": (1, 1_000_000, int),
        "min_depth": (1, 1_000_000, int),
        "min_variant_quality": (0, 10_000, float),
        "primer_min_length": (1, 1_000_000, int),
        "primer_min_quality": (0, 93, int),
    }
    for key, (lower, upper, converter) in numeric_rules.items():
        if key not in workflow:
            continue
        try:
            value = converter(workflow[key])
        except (TypeError, ValueError) as exc:
            raise ConfigError(f"workflow.{key} must be numeric") from exc
        if not lower <= value <= upper:
            raise ConfigError(f"workflow.{key} must be between {lower} and {upper}")

    minor = config.get("minor_variants", {})
    if not isinstance(minor, dict):
        raise ConfigError("minor_variants must be a YAML mapping")
    frequency = float(minor.get("min_frequency", 0.03))
    if not 0.001 <= frequency <= 0.5:
        raise ConfigError("minor_variants.min_frequency must be between 0.001 and 0.5")
    for key in ("min_depth", "min_base_quality"):
        if int(minor.get(key, 1)) < 1:
            raise ConfigError(f"minor_variants.{key} must be positive")

    profile_value = config.get("analysis", {}).get("profile")
    advisories: list[str] = []
    if profile_value:
        profile_path = resolve_path(config, profile_value)
        assert profile_path is not None
        if check_files and not profile_path.exists():
            raise ConfigError(f"analysis profile does not exist: {profile_path}")
        if profile_path.exists():
            validate_profile(load_yaml(profile_path), check_files=check_files)
    else:
        advisories.append("No analysis.profile configured; comparative antigenic report is disabled.")
    if not optional_files["nextclade_dataset"]:
        advisories.append("No Nextclade dataset configured; clade assignment is disabled.")
    return advisories


def validate_profile(profile: dict[str, Any], check_files: bool = True) -> None:
    reference = resolve_path(profile, profile.get("reference_fasta"))
    if reference is None:
        raise ConfigError("reference_fasta is required in the analysis profile")
    if check_files and not reference.exists():
        raise ConfigError(f"profile reference FASTA does not exist: {reference}")
    for key in ("vaccine_reference_fasta", "lineage_references"):
        path = resolve_path(profile, profile.get(key))
        if check_files and path and not path.exists():
            raise ConfigError(f"profile {key} does not exist: {path}")
    regions = parse_regions(profile)
    region_names = {region.name for region in regions}
    unknown_antigenic = set(profile.get("antigenic_sites", {})) - region_names
    unknown_markers = set(profile.get("escape_markers", {})) - region_names
    if unknown_antigenic or unknown_markers:
        unknown = sorted(unknown_antigenic | unknown_markers)
        raise ConfigError("profile references unknown protein(s): " + ", ".join(unknown))
    parse_markers(profile)
