from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from Bio import Align
from Bio.Seq import Seq

from .config import ConfigError, parse_markers, parse_regions, resolve_path
from .fasta import read_fasta
from .models import ProteinResult, SampleResult

VALID_BASES = set("ACGT")
MUTATION_PATTERN = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")


def _align(reference: str, query: str) -> dict[str, Any]:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    alignment = aligner.align(reference, query)[0]
    coordinates = alignment.coordinates

    projected = ["-"] * len(reference)
    changes: list[str] = []
    matches = 0
    callable_sites = 0
    insertions = 0

    for block in range(coordinates.shape[1] - 1):
        r0, r1 = int(coordinates[0, block]), int(coordinates[0, block + 1])
        q0, q1 = int(coordinates[1, block]), int(coordinates[1, block + 1])
        r_step, q_step = r1 - r0, q1 - q0
        if r_step and q_step:
            for offset in range(min(r_step, q_step)):
                r_index, q_index = r0 + offset, q0 + offset
                r_base, q_base = reference[r_index], query[q_index]
                projected[r_index] = q_base
                if r_base in VALID_BASES and q_base in VALID_BASES:
                    callable_sites += 1
                    if r_base == q_base:
                        matches += 1
                    else:
                        changes.append(f"{r_base}{r_index + 1}{q_base}")
        elif r_step:
            for r_index in range(r0, r1):
                if reference[r_index] in VALID_BASES:
                    callable_sites += 1
                changes.append(f"{reference[r_index]}{r_index + 1}-")
        elif q_step:
            inserted = query[q0:q1]
            insertions += len(inserted)
            changes.append(f"ins{r0}:{inserted}")

    return {
        "identity": round(matches / callable_sites, 6) if callable_sites else 0.0,
        "callable_sites": callable_sites,
        "reference_length": len(reference),
        "query_length": len(query),
        "n_content": round(query.count("N") / len(query), 6) if query else 1.0,
        "inserted_bases": insertions,
        "nucleotide_changes": changes,
        "projected_query": "".join(projected),
    }


def _translate_region(sequence: str, start: int, end: int, strand: int) -> str:
    nucleotides = sequence[start - 1 : end].replace("-", "N")
    if strand == -1:
        nucleotides = str(Seq(nucleotides).reverse_complement())
    nucleotides = nucleotides[: len(nucleotides) - (len(nucleotides) % 3)]
    return str(Seq(nucleotides).translate())


def _amino_acid_changes(reference: str, query: str) -> list[str]:
    changes: list[str] = []
    for position, (ref_aa, query_aa) in enumerate(zip(reference, query), start=1):
        if query_aa == "X":
            continue
        if ref_aa != query_aa:
            changes.append(f"{ref_aa}{position}{query_aa}")
    if len(query) < len(reference):
        changes.extend(f"{aa}{pos}-" for pos, aa in enumerate(reference[len(query) :], len(query) + 1))
    return changes


def _load_lineages(profile: dict[str, Any]) -> dict[str, dict[str, str]]:
    path = resolve_path(profile, profile.get("lineage_references"))
    if path is None:
        return {}
    records = read_fasta(path)
    lineages: dict[str, dict[str, str]] = {}
    for header, sequence in records.items():
        try:
            lineage, segment = header.split("|", 1)
        except ValueError as exc:
            raise ConfigError(
                "lineage reference FASTA IDs must use the form lineage|segment"
            ) from exc
        lineages.setdefault(lineage, {})[segment] = sequence
    return lineages


def _genotype_and_shift(
    queries: dict[str, str],
    lineages: dict[str, dict[str, str]],
    segmented: bool,
    min_margin: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    segment_calls: dict[str, dict[str, Any]] = {}
    for segment, query in queries.items():
        scores = []
        for lineage, sequences in lineages.items():
            if segment in sequences:
                scores.append((lineage, _align(sequences[segment], query)["identity"]))
        scores.sort(key=lambda item: item[1], reverse=True)
        if not scores:
            segment_calls[segment] = {"lineage": None, "identity": None, "margin": None}
            continue
        best_lineage, best_score = scores[0]
        runner_up = scores[1][1] if len(scores) > 1 else 0.0
        margin = best_score - runner_up
        segment_calls[segment] = {
            "lineage": best_lineage,
            "identity": round(best_score, 6),
            "margin": round(margin, 6),
            "confident": margin >= min_margin,
        }

    confident = {
        segment: call["lineage"]
        for segment, call in segment_calls.items()
        if call.get("confident") and call.get("lineage")
    }
    lineages_used = sorted(set(confident.values()))
    candidate = segmented and len(lineages_used) >= 2
    genotype = {
        "method": "best-reference identity",
        "segment_calls": segment_calls,
        "composite": "/".join(
            f"{segment}:{call['lineage'] or 'unassigned'}"
            for segment, call in sorted(segment_calls.items())
        ),
    }
    shift = {
        "applicable": segmented,
        "candidate_reassortment_signal": candidate,
        "confident_lineages": lineages_used,
        "interpretation": (
            "Segments have confident best matches to multiple reference lineages; "
            "confirm with phylogenetics, contamination checks, and epidemiologic context."
            if candidate
            else "No multi-lineage segment pattern passed the configured identity-margin screen."
            if segmented
            else "Not evaluated: the profile is not configured as a segmented virus."
        ),
    }
    return genotype, shift


def analyze_genome(
    sample: str,
    query_fasta: Path,
    profile: dict[str, Any],
) -> SampleResult:
    reference_path = resolve_path(profile, profile.get("reference_fasta"))
    if reference_path is None:
        raise ConfigError("reference_fasta is required in the analysis profile")
    references = read_fasta(reference_path)
    queries = read_fasta(query_fasta)
    missing = sorted(set(references) - set(queries))
    if missing:
        raise ConfigError(f"query FASTA is missing reference segment(s): {', '.join(missing)}")

    segment_results: dict[str, dict[str, Any]] = {}
    projected: dict[str, str] = {}
    for segment, reference in references.items():
        result = _align(reference, queries[segment])
        projected[segment] = result.pop("projected_query")
        segment_results[segment] = result

    regions = parse_regions(profile)
    antigenic_sites = {
        str(name): {int(position) for position in positions}
        for name, positions in profile.get("antigenic_sites", {}).items()
    }
    markers = parse_markers(profile)
    protein_results: dict[str, ProteinResult] = {}
    total_antigenic_sites = 0
    changed_antigenic_sites = 0
    total_aa = 0
    changed_aa = 0
    total_marker_weight = sum(m.weight for values in markers.values() for m in values)
    matched_marker_weight = 0.0

    for region in regions:
        if region.segment not in references:
            raise ConfigError(f"coding region {region.name} names unknown segment {region.segment}")
        if region.end > len(references[region.segment]):
            raise ConfigError(f"coding region {region.name} extends beyond its reference segment")
        ref_protein = _translate_region(
            references[region.segment], region.start, region.end, region.strand
        )
        query_protein = _translate_region(
            projected[region.segment], region.start, region.end, region.strand
        )
        aa_changes = _amino_acid_changes(ref_protein, query_protein)
        change_positions = {
            int(match.group(2))
            for change in aa_changes
            if (match := MUTATION_PATTERN.match(change))
        }
        sites = antigenic_sites.get(region.name, set())
        site_changes = [
            change
            for change in aa_changes
            if (match := MUTATION_PATTERN.match(change)) and int(match.group(2)) in sites
        ]
        matched = []
        for marker in markers.get(region.name, []):
            if marker.mutation in aa_changes:
                matched.append(
                    {
                        "mutation": marker.mutation,
                        "weight": marker.weight,
                        "evidence": marker.evidence,
                    }
                )
                matched_marker_weight += marker.weight
        protein_results[region.name] = ProteinResult(
            name=region.name,
            amino_acid_changes=aa_changes,
            antigenic_site_changes=site_changes,
            matched_markers=matched,
        )
        total_antigenic_sites += len(sites)
        changed_antigenic_sites += len(change_positions & sites)
        total_aa += len(ref_protein)
        changed_aa += len(aa_changes)

    antigenic_rate = changed_antigenic_sites / total_antigenic_sites if total_antigenic_sites else 0
    aa_divergence = changed_aa / total_aa if total_aa else 0
    drift_index = min(100.0, 60 * antigenic_rate + 40 * min(1.0, aa_divergence / 0.05))
    drift = {
        "evidence_index": round(drift_index, 2),
        "changed_antigenic_sites": changed_antigenic_sites,
        "configured_antigenic_sites": total_antigenic_sites,
        "amino_acid_divergence": round(aa_divergence, 6),
        "formula": "60*antigenic-site change fraction + 40*min(1, AA divergence/0.05)",
    }

    vaccine_path = resolve_path(profile, profile.get("vaccine_reference_fasta"))
    vaccine_identity = None
    if vaccine_path:
        vaccine_records = read_fasta(vaccine_path)
        identities = [
            _align(vaccine_records[segment], queries[segment])["identity"]
            for segment in references
            if segment in vaccine_records
        ]
        if identities:
            vaccine_identity = sum(identities) / len(identities)

    marker_component = (
        70 * matched_marker_weight / total_marker_weight if total_marker_weight else 0
    )
    antigenic_component = 20 * antigenic_rate
    distance_component = (
        10 * min(1.0, (1 - vaccine_identity) / 0.05) if vaccine_identity is not None else 0
    )
    escape_score = min(100.0, marker_component + antigenic_component + distance_component)
    priority = "high" if escape_score >= 50 else "moderate" if escape_score >= 25 else "low"
    escape = {
        "priority_score": round(escape_score, 2),
        "priority_band": priority,
        "matched_marker_weight": round(matched_marker_weight, 3),
        "configured_marker_weight": round(total_marker_weight, 3),
        "vaccine_nucleotide_identity": (
            round(vaccine_identity, 6) if vaccine_identity is not None else None
        ),
        "components": {
            "marker_evidence_max_70": round(marker_component, 2),
            "antigenic_sites_max_20": round(antigenic_component, 2),
            "vaccine_distance_max_10": round(distance_component, 2),
        },
        "interpretation": (
            "Heuristic triage score only; it is not a neutralization estimate, vaccine-effectiveness "
            "prediction, or clinical conclusion."
        ),
    }

    lineages = _load_lineages(profile)
    genotype, shift = _genotype_and_shift(
        queries,
        lineages,
        bool(profile.get("segmented", False)),
        float(profile.get("lineage_min_margin", 0.002)),
    )
    warnings = []
    if any(value["n_content"] > 0.05 for value in segment_results.values()):
        warnings.append("At least one segment has >5% N bases; interpret calls cautiously.")
    if not lineages:
        warnings.append("No lineage_references configured; genotype calls are unassigned.")
    if not total_antigenic_sites:
        warnings.append("No antigenic sites configured; drift evidence is sequence-divergence only.")
    if not total_marker_weight:
        warnings.append("No evidence-linked escape markers configured.")

    return SampleResult(
        sample=sample,
        segments=segment_results,
        proteins=protein_results,
        genotype=genotype,
        drift=drift,
        shift_screen=shift,
        vaccine_escape_screen=escape,
        warnings=warnings,
    )
