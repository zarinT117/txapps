from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path
from typing import Any

from . import __author__, __author_email__, __version__
from .analysis import analyze_genome
from .config import ConfigError, load_samples, load_yaml, resolve_path, validate_config
from .provenance import collect_tool_versions, sha256_file
from .qc import write_qc_summary
from .report import write_reports
from .runner import Runner

CORE_TOOLS = ("minimap2", "samtools", "bcftools")
PLATFORM_TOOLS = {
    "illumina": ("fastp", "megahit"),
    "nanopore": ("filtlong", "flye"),
}
OPTIONAL_TOOLS = ("ivar", "nextclade")


def tool_status(include_nextclade: bool = True) -> dict[str, str | None]:
    tools = list(CORE_TOOLS)
    tools.extend(tool for values in PLATFORM_TOOLS.values() for tool in values)
    tools.append("ivar")
    if include_nextclade:
        tools.append("nextclade")
    return {tool: shutil.which(tool) for tool in dict.fromkeys(tools)}


def required_tools_for_config(config: dict[str, Any]) -> list[str]:
    samples = load_samples(config)
    tools = list(CORE_TOOLS)
    for platform in sorted({sample.platform for sample in samples}):
        tools.extend(PLATFORM_TOOLS[platform])
    workflow = config.get("workflow", {})
    if workflow.get("primer_bed") or config.get("minor_variants", {}).get("enabled"):
        tools.append("ivar")
    if workflow.get("nextclade_dataset"):
        tools.append("nextclade")
    return list(dict.fromkeys(tools))


def _low_coverage_bed(depth_path: Path, bed_path: Path, min_depth: int) -> None:
    intervals: list[tuple[str, int, int]] = []
    active_contig: str | None = None
    active_start: int | None = None
    previous_pos = 0
    with depth_path.open(encoding="utf-8") as handle:
        for line in handle:
            contig, pos_text, depth_text = line.rstrip().split("\t")[:3]
            position, depth = int(pos_text), int(depth_text)
            is_low = depth < min_depth
            if is_low and (active_contig != contig or active_start is None):
                if active_contig is not None and active_start is not None:
                    intervals.append((active_contig, active_start - 1, previous_pos))
                active_contig, active_start = contig, position
            elif not is_low and active_contig is not None and active_start is not None:
                intervals.append((active_contig, active_start - 1, previous_pos))
                active_contig, active_start = None, None
            previous_pos = position
    if active_contig is not None and active_start is not None:
        intervals.append((active_contig, active_start - 1, previous_pos))
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    with bed_path.open("w", encoding="utf-8", newline="\n") as handle:
        for contig, start, end in intervals:
            handle.write(f"{contig}\t{start}\t{end}\n")


def _existing_hashes(paths: list[Path | None]) -> dict[str, str]:
    return {str(path): sha256_file(path) for path in paths if path and path.exists()}


def _deplete_host(
    runner: Runner,
    sample_name: str,
    platform: str,
    reads: list[Path],
    host_reference: Path,
    qc_dir: Path,
    threads: int,
) -> list[Path]:
    depleted_dir = qc_dir / "host_depleted"
    depleted_dir.mkdir(parents=True, exist_ok=True)
    preset = "sr" if platform == "illumina" else "map-ont"
    host_bam = depleted_dir / "host_screen.bam"
    runner.pipe(
        [
            ["minimap2", "-t", str(threads), "-ax", preset, host_reference, *reads],
            ["samtools", "sort", "-@", str(threads), "-o", host_bam],
        ],
        "02a_host_screen",
        expected_outputs=[host_bam],
    )

    if platform == "illumina" and len(reads) == 2:
        unmapped_bam = depleted_dir / "unmapped.namesort.bam"
        depleted_r1 = depleted_dir / f"{sample_name}_host_depleted_R1.fastq.gz"
        depleted_r2 = depleted_dir / f"{sample_name}_host_depleted_R2.fastq.gz"
        runner.pipe(
            [
                ["samtools", "view", "-@", str(threads), "-b", "-f", "12", "-F", "256", host_bam],
                ["samtools", "sort", "-@", str(threads), "-n", "-o", unmapped_bam],
            ],
            "02b_host_unmapped_pairs",
            expected_outputs=[unmapped_bam],
        )
        runner.run(
            [
                "samtools",
                "fastq",
                "-@",
                str(threads),
                "-n",
                "-1",
                depleted_r1,
                "-2",
                depleted_r2,
                "-0",
                os.devnull,
                "-s",
                os.devnull,
                unmapped_bam,
            ],
            "02c_host_depleted_fastq",
            expected_outputs=[depleted_r1, depleted_r2],
        )
        return [depleted_r1, depleted_r2]

    unmapped_bam = depleted_dir / "unmapped.bam"
    depleted_reads = depleted_dir / f"{sample_name}_host_depleted.fastq.gz"
    runner.run(
        ["samtools", "view", "-@", str(threads), "-b", "-f", "4", "-F", "256", "-o", unmapped_bam, host_bam],
        "02b_host_unmapped_reads",
        expected_outputs=[unmapped_bam],
    )
    runner.run(
        ["samtools", "fastq", "-@", str(threads), "-n", "-0", depleted_reads, "-s", os.devnull, unmapped_bam],
        "02c_host_depleted_fastq",
        expected_outputs=[depleted_reads],
    )
    return [depleted_reads]


def _trim_primers(
    runner: Runner,
    bam: Path,
    primer_bed: Path,
    mapping_dir: Path,
    min_length: int,
    min_quality: int,
    threads: int,
) -> Path:
    prefix = mapping_dir / "primer_trimmed"
    trimmed_bam = prefix.with_suffix(".bam")
    sorted_bam = mapping_dir / "primer_trimmed.sorted.bam"
    runner.run(
        ["ivar", "trim", "-i", bam, "-b", primer_bed, "-p", prefix, "-m", str(min_length), "-q", str(min_quality)],
        "04a_trim_primers",
        expected_outputs=[trimmed_bam],
    )
    runner.run(
        ["samtools", "sort", "-@", str(threads), "-o", sorted_bam, trimmed_bam],
        "04b_sort_primer_trimmed_bam",
        expected_outputs=[sorted_bam],
    )
    runner.run(["samtools", "index", "-@", str(threads), sorted_bam], "04c_index_primer_trimmed_bam")
    return sorted_bam


def _write_cohort_summary(sample_dirs: list[Path], output_path: Path) -> None:
    fields = [
        "sample",
        "result_directory",
        "mean_depth",
        "breadth_at_10x",
        "n_content_max",
        "drift_index",
        "escape_priority",
        "escape_score",
        "candidate_reassortment_signal",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for sample_dir in sample_dirs:
            qc_path = sample_dir / "05_report" / "qc_summary.json"
            analysis_path = sample_dir / "05_report" / "analysis.json"
            qc = json.loads(qc_path.read_text(encoding="utf-8")) if qc_path.exists() else {}
            analysis = (
                json.loads(analysis_path.read_text(encoding="utf-8")) if analysis_path.exists() else {}
            )
            segments = analysis.get("segments", {})
            n_content_max = max((values.get("n_content", 0) for values in segments.values()), default="")
            writer.writerow(
                {
                    "sample": sample_dir.name,
                    "result_directory": sample_dir,
                    "mean_depth": qc.get("coverage", {}).get("mean_depth", ""),
                    "breadth_at_10x": qc.get("coverage", {}).get("breadth_at_10x", ""),
                    "n_content_max": n_content_max,
                    "drift_index": analysis.get("drift", {}).get("evidence_index", ""),
                    "escape_priority": analysis.get("vaccine_escape_screen", {}).get("priority_band", ""),
                    "escape_score": analysis.get("vaccine_escape_screen", {}).get("priority_score", ""),
                    "candidate_reassortment_signal": analysis.get("shift_screen", {}).get(
                        "candidate_reassortment_signal", ""
                    ),
                }
            )


def run_pipeline(config_path: Path, dry_run: bool = False, resume: bool = False) -> list[Path]:
    config = load_yaml(config_path)
    advisories = validate_config(config, check_files=not dry_run)
    for advisory in advisories:
        print(f"[advisory] {advisory}")

    required_tools = required_tools_for_config(config)
    if not dry_run:
        missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing_tools:
            raise ConfigError(
                "required external tool(s) missing: "
                + ", ".join(missing_tools)
                + ". Install them with the documented conda environment."
            )

    workflow = config.get("workflow", {})
    samples = load_samples(config)
    reference = resolve_path(config, workflow.get("reference_fasta"))
    if reference is None:
        raise ConfigError("workflow.reference_fasta is required")
    host_reference = resolve_path(config, workflow.get("host_reference_fasta"))
    primer_bed = resolve_path(config, workflow.get("primer_bed"))
    nextclade_dataset = resolve_path(config, workflow.get("nextclade_dataset"))
    profile_path = resolve_path(config, config.get("analysis", {}).get("profile"))
    output_root = resolve_path(config, workflow.get("output_dir", "results"))
    assert output_root is not None

    threads = int(workflow.get("threads", 4))
    min_depth = int(workflow.get("min_depth", 10))
    min_quality = float(workflow.get("min_variant_quality", 20))
    min_read_length = int(workflow.get("min_read_length", 50))
    primer_min_length = int(workflow.get("primer_min_length", min_read_length))
    primer_min_quality = int(workflow.get("primer_min_quality", 20))
    minor = config.get("minor_variants", {})
    minor_enabled = bool(minor.get("enabled", False))
    minor_min_frequency = float(minor.get("min_frequency", 0.03))
    minor_min_depth = int(minor.get("min_depth", min_depth))
    minor_min_base_quality = int(minor.get("min_base_quality", 20))
    annotation_gff = resolve_path(config, workflow.get("annotation_gff"))

    completed: list[Path] = []
    if not dry_run:
        first_log_dir = output_root / "_logs"
        Runner(first_log_dir, dry_run=False, resume=resume).run(
            ["samtools", "faidx", reference],
            "00_index_reference",
            expected_outputs=[Path(str(reference) + ".fai")],
        )

    for sample in samples:
        sample_dir = output_root / sample.name
        qc_dir = sample_dir / "01_qc"
        assembly_dir = sample_dir / "02_assembly"
        mapping_dir = sample_dir / "03_mapping"
        genotype_dir = sample_dir / "04_genotype"
        report_dir = sample_dir / "05_report"
        log_dir = sample_dir / "logs"
        for directory in (qc_dir, assembly_dir, mapping_dir, genotype_dir, report_dir, log_dir):
            directory.mkdir(parents=True, exist_ok=True)
        runner = Runner(log_dir, dry_run=dry_run, resume=resume)

        if sample.platform == "illumina":
            clean_r1 = qc_dir / "clean_R1.fastq.gz"
            clean_r2 = qc_dir / "clean_R2.fastq.gz" if sample.r2 else None
            fastp_command: list[str | Path] = [
                "fastp",
                "--in1",
                sample.r1,
                "--out1",
                clean_r1,
                "--json",
                qc_dir / "fastp.json",
                "--html",
                qc_dir / "fastp.html",
                "--thread",
                str(threads),
                "--length_required",
                str(min_read_length),
            ]
            if sample.r2 and clean_r2:
                fastp_command += ["--in2", sample.r2, "--out2", clean_r2, "--detect_adapter_for_pe"]
            runner.run(fastp_command, "01_fastp", expected_outputs=[clean_r1] + ([clean_r2] if clean_r2 else []))
            reads = [clean_r1] + ([clean_r2] if clean_r2 else [])
            preset = "sr"
        else:
            clean_r1 = qc_dir / "filtered.fastq.gz"
            runner.run(
                ["filtlong", "--min_length", str(min_read_length), sample.r1],
                "01_filtlong",
                stdout_path=clean_r1,
                expected_outputs=[clean_r1],
            )
            reads = [clean_r1]
            preset = "map-ont"

        if host_reference:
            reads = _deplete_host(runner, sample.name, sample.platform, reads, host_reference, qc_dir, threads)

        if sample.platform == "illumina":
            runner.run(
                ["megahit", "--out-dir", assembly_dir, "--num-cpu-threads", str(threads)]
                + (["-1", reads[0], "-2", reads[1]] if len(reads) == 2 else ["-r", reads[0]]),
                "02_assemble_megahit",
                expected_outputs=[assembly_dir / "final.contigs.fa"],
            )
            contigs = assembly_dir / "final.contigs.fa"
        else:
            runner.run(
                ["flye", "--nano-hq", reads[0], "--out-dir", assembly_dir, "--threads", str(threads)],
                "02_assemble_flye",
                expected_outputs=[assembly_dir / "assembly.fasta"],
            )
            contigs = assembly_dir / "assembly.fasta"

        bam = mapping_dir / "aligned.bam"
        runner.pipe(
            [
                ["minimap2", "-t", str(threads), "-ax", preset, reference, *reads],
                ["samtools", "sort", "-@", str(threads), "-o", bam],
            ],
            "03_map_and_sort",
            expected_outputs=[bam],
        )
        runner.run(["samtools", "index", "-@", str(threads), bam], "04_index_bam")
        runner.run(["samtools", "flagstat", "-@", str(threads), bam], "04_flagstat", stdout_path=mapping_dir / "flagstat.txt")

        call_bam = (
            _trim_primers(runner, bam, primer_bed, mapping_dir, primer_min_length, primer_min_quality, threads)
            if primer_bed
            else bam
        )

        raw_vcf = mapping_dir / "variants.raw.vcf.gz"
        normalized_vcf = mapping_dir / "variants.normalized.vcf.gz"
        filtered_vcf = mapping_dir / "variants.filtered.vcf.gz"
        runner.pipe(
            [
                [
                    "bcftools",
                    "mpileup",
                    "--threads",
                    str(threads),
                    "-Ou",
                    "-f",
                    reference,
                    "-a",
                    "INFO/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP",
                    call_bam,
                ],
                ["bcftools", "call", "--threads", str(threads), "--ploidy", "1", "-mv", "-Oz", "-o", raw_vcf],
            ],
            "05_call_variants",
            expected_outputs=[raw_vcf],
        )
        runner.run(
            ["bcftools", "norm", "--threads", str(threads), "-f", reference, "-m", "-any", "-Oz", "-o", normalized_vcf, raw_vcf],
            "06_normalize_variants",
            expected_outputs=[normalized_vcf],
        )
        runner.run(
            [
                "bcftools",
                "filter",
                "--threads",
                str(threads),
                "-i",
                f"QUAL>={min_quality} && INFO/DP>={min_depth}",
                "-Oz",
                "-o",
                filtered_vcf,
                normalized_vcf,
            ],
            "07_filter_variants",
            expected_outputs=[filtered_vcf],
        )
        runner.run(["bcftools", "index", "-f", filtered_vcf], "08_index_vcf")
        runner.run(["bcftools", "stats", filtered_vcf], "08_vcf_stats", stdout_path=mapping_dir / "variants.stats.txt")

        depth_path = mapping_dir / "depth.tsv"
        runner.run(["samtools", "depth", "-aa", call_bam], "09_depth", stdout_path=depth_path, expected_outputs=[depth_path])
        mask_bed = mapping_dir / "low_coverage.bed"
        if not dry_run:
            _low_coverage_bed(depth_path, mask_bed, min_depth)
        consensus = mapping_dir / "consensus.fasta"
        runner.run(
            ["bcftools", "consensus", "-f", reference, "--mask", mask_bed, "--missing", "N", filtered_vcf],
            "10_consensus",
            stdout_path=consensus,
            expected_outputs=[consensus],
        )

        if minor_enabled:
            minor_prefix = genotype_dir / "minor_variants"
            ivar_command: list[str | Path] = [
                "ivar",
                "variants",
                "-p",
                minor_prefix,
                "-q",
                str(minor_min_base_quality),
                "-t",
                str(minor_min_frequency),
                "-m",
                str(minor_min_depth),
                "-r",
                reference,
            ]
            if annotation_gff:
                ivar_command.extend(["-g", annotation_gff])
            runner.pipe(
                [["samtools", "mpileup", "-aa", "-A", "-d", "0", "-B", "-Q", "0", "-f", reference, call_bam], ivar_command],
                "11_minor_variants",
                expected_outputs=[minor_prefix.with_suffix(".tsv")],
            )

        if nextclade_dataset:
            runner.run(
                ["nextclade", "run", "--input-dataset", nextclade_dataset, "--output-all", genotype_dir, consensus],
                "12_nextclade",
                expected_outputs=[genotype_dir / "nextclade.tsv"],
            )

        if not dry_run:
            qc_payload = write_qc_summary(
                report_dir / "qc_summary.json",
                contigs,
                depth_path,
                additional={
                    "reference": str(reference),
                    "host_depletion": bool(host_reference),
                    "primer_trimming": bool(primer_bed),
                    "consensus_fasta": str(consensus),
                    "filtered_vcf": str(filtered_vcf),
                },
            )
            if profile_path:
                profile = load_yaml(profile_path)
                result = analyze_genome(sample.name, consensus, profile)
                write_reports(result, report_dir)
            manifest_payload = {
                "software": {"name": "ViroFlow", "version": __version__, "author": __author__, "author_email": __author_email__},
                "sample": sample.name,
                "platform": sample.platform,
                "inputs": _existing_hashes([sample.r1, sample.r2, reference, host_reference, primer_bed, profile_path]),
                "outputs": {
                    "contigs": str(contigs),
                    "bam": str(call_bam),
                    "filtered_vcf": str(filtered_vcf),
                    "consensus": str(consensus),
                    "qc_summary": str(report_dir / "qc_summary.json"),
                },
                "parameters": {
                    "min_depth": min_depth,
                    "min_variant_quality": min_quality,
                    "min_read_length": min_read_length,
                    "minor_variants": minor if minor_enabled else {"enabled": False},
                },
                "tool_versions": collect_tool_versions(required_tools),
                "qc": qc_payload,
            }
        else:
            manifest_payload = {
                "software": {"name": "ViroFlow", "version": __version__, "author": __author__, "author_email": __author_email__},
                "sample": sample.name,
                "platform": sample.platform,
                "reference": str(reference),
                "dry_run_only": True,
            }

        runner.write_manifest(sample_dir / "run_manifest.json", manifest_payload)
        completed.append(sample_dir)

    if completed and not dry_run:
        _write_cohort_summary(completed, completed[0].parent / "cohort_summary.tsv")
    return completed


def write_summary(sample_dirs: list[Path], output_path: Path) -> None:
    _write_cohort_summary(sample_dirs, output_path)
