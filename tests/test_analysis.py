from pathlib import Path

import pytest

from viroflow.analysis import _align, analyze_genome
from viroflow.config import load_yaml

DEMO = Path(__file__).parent / "fixtures" / "synthetic"


def test_alignment_reports_substitution_and_missing_data():
    result = _align("ACGT", "ATNT")
    assert result["nucleotide_changes"] == ["C2T"]
    assert result["callable_sites"] == 3
    assert result["identity"] == pytest.approx(2 / 3, abs=1e-6)


def test_demo_detects_mutations_and_reassortment_candidate():
    profile = load_yaml(DEMO / "profile.yaml")
    result = analyze_genome("demo", DEMO / "query.fasta", profile)
    assert "K2N" in result.proteins["HA"].amino_acid_changes
    assert "F4Y" in result.proteins["HA"].antigenic_site_changes
    assert result.shift_screen["candidate_reassortment_signal"] is True
    assert result.genotype["segment_calls"]["HA"]["lineage"] == "lineage_A"
    assert result.genotype["segment_calls"]["NA"]["lineage"] == "lineage_B"
    assert result.vaccine_escape_screen["priority_score"] > 0


def test_identical_sequences_have_no_changes(tmp_path):
    profile = load_yaml(DEMO / "profile.yaml")
    result = analyze_genome("reference", DEMO / "reference.fasta", profile)
    assert not result.segments["HA"]["nucleotide_changes"]
    assert not result.proteins["HA"].amino_acid_changes
