from pathlib import Path

from viroflow.pipeline import _low_coverage_bed
from viroflow.runner import display_command


def test_low_coverage_positions_are_merged_into_bed(tmp_path: Path):
    depth = tmp_path / "depth.tsv"
    depth.write_text("seg\t1\t0\nseg\t2\t3\nseg\t3\t12\nseg\t4\t2\n", encoding="utf-8")
    bed = tmp_path / "mask.bed"
    _low_coverage_bed(depth, bed, min_depth=10)
    assert bed.read_text(encoding="utf-8") == "seg\t0\t2\nseg\t3\t4\n"


def test_commands_are_shell_escaped_for_manifest_display():
    rendered = display_command(["bcftools", "filter", "-i", "QUAL>=20 && INFO/DP>=10"])
    assert "bcftools filter" in rendered
    assert "QUAL>=20 && INFO/DP>=10" in rendered

