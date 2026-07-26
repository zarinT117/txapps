from __future__ import annotations

import csv
import shutil
from pathlib import Path

from .analysis import analyze_genome
from .config import ConfigError, load_samples, load_yaml, resolve_path
from .report import write_reports
from .runner import Runner

REQUIRED_TOOLS = ("fastp", "megahit", "filtlong", "flye", "minimap2", "samtools", "bcftools")


def tool_status(include_nextclade: bool = True) -> dict[str, str | None]:
    tools = list(REQUIRED_TOOLS)
    if include_nextclade:
        tools.append("nextclade")
    return {tool: shutil.which(tool) for tool in tools}


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
    with bed_path.open("w", encoding="utf-8", newline="\n") as handle:
        for contig, start, end in intervals:
            handle.write(f"{contig}\t{start}\t{end}\n")


def run_pipeline(config_path: Path, dry_run: bool = False) -> list[Path]:
    config = load_yaml(config_path)
    workflow = config.get("workflow", {})
    samples = load_samples(config)
    reference = resolve_path(config, workflow.get("reference_fasta"))
    if reference is None:
        raise ConfigError("workflow.reference_fasta is required")
    if not dry_run:
        if not reference.exists():
            raise ConfigError(f"reference FASTA does not exist: {reference}")
        missing_inputs = [
            str(path)
            for sample in samples
            for path in (sample.r1, sample.r2)
            if path is not None and not path.exists()
        ]
        if missing_inputs:
            raise ConfigError("read file(s) do not exist: " + ", ".join(missing_inputs))

    output_root = resolve_path(config, workflow.get("output_dir", "results"))
    assert output_root is not None
    threads = int(workflow.get("threads", 4))
    min_depth = int(workflow.get("min_depth", 10))
    min_quality = float(workflow.get("min_variant_quality", 20))
    min_read_length = int(workflow.get("min_read_length", 50))
    nextclade_dataset = resolve_path(config, workflow.get("nextclade_dataset"))
    profile_path = resolve_path(config, config.get("analysis", {}).get("profile"))
    completed: list[Path] = []

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
        runner = Runner(log_dir, dry_run=dry_run)

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
            runner.run(fastp_command, "01_fastp")
            megahit_command: list[str | Path] = [
                "megahit",
                "--out-dir",
                assembly_dir,
                "--num-cpu-threads",
                str(threads),
            ]
            if clean_r2:
                megahit_command += ["-1", clean_r1, "-2", clean_r2]
            else:
                megahit_command += ["-r", clean_r1]
            runner.run(megahit_command, "02_megahit")
            contigs = assembly_dir / "final.contigs.fa"
            reads = [clean_r1] + ([clean_r2] if clean_r2 else [])
            preset = "sr"
        else:
            clean_r1 = qc_dir / "filtered.fastq"
            runner.run(
                ["filtlong", "--min_length", str(min_read_length), sample.r1],
                "01_filtlong",
                stdout_path=clean_r1,
            )
            runner.run(
                [
                    "flye",
                    "--nano-hq",
                    clean_r1,
                    "--out-dir",
                    assembly_dir,
                    "--threads",
                    str(threads),
                ],
                "02_flye",
            )
            contigs = assembly_dir / "assembly.fasta"
            reads = [clean_r1]
            preset = "map-ont"

        bam = mapping_dir / "aligned.bam"
        runner.pipe(
            [
                ["minimap2", "-t", str(threads), "-ax", preset, reference, *reads],
                ["samtools", "sort", "-@", str(threads), "-o", bam],
            ],
            "03_map_and_sort",
        )
        runner.run(["samtools", "index", "-@", str(threads), bam], "04_index_bam")

        raw_vcf = mapping_dir / "variants.raw.vcf.gz"
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
                    "INFO/DP,FORMAT/DP",
                    bam,
                ],
                [
                    "bcftools",
                    "call",
                    "--threads",
                    str(threads),
                    "--ploidy",
                    "1",
                    "-mv",
                    "-Oz",
                    "-o",
                    raw_vcf,
                ],
            ],
            "05_call_variants",
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
                raw_vcf,
            ],
            "06_filter_variants",
        )
        runner.run(["bcftools", "index", "-f", filtered_vcf], "07_index_vcf")

        depth_path = mapping_dir / "depth.tsv"
        runner.run(["samtools", "depth", "-aa", bam], "08_depth", stdout_path=depth_path)
        mask_bed = mapping_dir / "low_coverage.bed"
        if not dry_run:
            _low_coverage_bed(depth_path, mask_bed, min_depth)
        consensus = mapping_dir / "consensus.fasta"
        runner.run(
            [
                "bcftools",
                "consensus",
                "-f",
                reference,
                "--mask",
                mask_bed,
                "--missing",
                "N",
                filtered_vcf,
            ],
            "09_consensus",
            stdout_path=consensus,
        )

        if nextclade_dataset:
            runner.run(
                [
                    "nextclade",
                    "run",
                    "--input-dataset",
                    nextclade_dataset,
                    "--output-all",
                    genotype_dir,
                    consensus,
                ],
                "10_nextclade",
            )

        if profile_path and not dry_run:
            profile = load_yaml(profile_path)
            result = analyze_genome(sample.name, consensus, profile)
            write_reports(result, report_dir)

        runner.write_manifest(
            sample_dir / "run_manifest.json",
            {
                "sample": sample.name,
                "platform": sample.platform,
                "reference": str(reference),
                "contigs": str(contigs),
                "consensus": str(consensus),
                "min_depth": min_depth,
                "min_variant_quality": min_quality,
            },
        )
        completed.append(sample_dir)
    return completed


def write_summary(sample_dirs: list[Path], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "result_directory"])
        for sample_dir in sample_dirs:
            writer.writerow([sample_dir.name, sample_dir])

