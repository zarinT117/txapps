from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Sample:
    name: str
    r1: Path
    r2: Path | None = None
    platform: str = "illumina"


@dataclass(frozen=True)
class CodingRegion:
    name: str
    segment: str
    start: int
    end: int
    strand: int = 1


@dataclass(frozen=True)
class Marker:
    mutation: str
    weight: float
    evidence: str


@dataclass
class ProteinResult:
    name: str
    amino_acid_changes: list[str] = field(default_factory=list)
    antigenic_site_changes: list[str] = field(default_factory=list)
    matched_markers: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class SampleResult:
    sample: str
    segments: dict[str, dict[str, Any]]
    proteins: dict[str, ProteinResult]
    genotype: dict[str, Any]
    drift: dict[str, Any]
    shift_screen: dict[str, Any]
    vaccine_escape_screen: dict[str, Any]
    warnings: list[str] = field(default_factory=list)

