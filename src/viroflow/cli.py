from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from . import __author__, __author_email__, __version__
from .analysis import analyze_genome
from .config import ConfigError, load_yaml, validate_config
from .ml import MLError, extract_features, predict_classifier, score_anomalies, train_classifier
from .pipeline import required_tools_for_config, run_pipeline, tool_status, write_summary
from .report import write_reports
from .templates import (
    CONFIG_TEMPLATE,
    PROFILE_TEMPLATE,
    SAMPLES_TEMPLATE,
    SARS_COV_2_CONFIG_TEMPLATE,
)


def banner() -> str:
    return f"ViroFlow {__version__} | Author: {__author__} <{__author_email__}>"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viroflow",
        description="Viral raw-read workflow and comparative antigenic evidence screen.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__} - author: {__author__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="write starter configuration files")
    init_parser.add_argument("directory", nargs="?", type=Path, default=Path("viroflow-project"))
    init_parser.add_argument("--preset", choices=("generic", "sars-cov-2"), default="generic")

    validate_parser = subparsers.add_parser("validate", help="validate a workflow configuration")
    validate_parser.add_argument("--config", type=Path, required=True)
    validate_parser.add_argument("--no-file-check", action="store_true")

    run_parser = subparsers.add_parser("run", help="run raw-read workflow")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--dry-run", action="store_true", help="print commands without executing")
    run_parser.add_argument("--resume", action="store_true", help="skip steps whose expected outputs already exist")

    analyze_parser = subparsers.add_parser(
        "analyze", help="analyze an assembled genome FASTA with an evidence profile"
    )
    analyze_parser.add_argument("--profile", type=Path, required=True)
    analyze_parser.add_argument("--input", type=Path, required=True)
    analyze_parser.add_argument("--sample", default="sample")
    analyze_parser.add_argument("--output", type=Path, required=True)

    doctor_parser = subparsers.add_parser("doctor", help="check external bioinformatics tools")
    doctor_parser.add_argument("--json", action="store_true", dest="as_json")
    doctor_parser.add_argument("--config", type=Path, help="check only tools required by this config")

    ml_parser = subparsers.add_parser("ml", help="machine-learning utilities for cohort prioritization")
    ml_sub = ml_parser.add_subparsers(dest="ml_command", required=True)

    features_parser = ml_sub.add_parser("features", help="extract ML-ready features from ViroFlow results")
    features_parser.add_argument("--results", type=Path, required=True)
    features_parser.add_argument("--output", type=Path, required=True)

    train_parser = ml_sub.add_parser("train", help="train a labelled escape/drift classifier")
    train_parser.add_argument("--features", type=Path, required=True)
    train_parser.add_argument("--label-column", required=True)
    train_parser.add_argument("--model", type=Path, required=True)
    train_parser.add_argument("--report", type=Path, required=True)
    train_parser.add_argument("--sample-column", default="sample")

    predict_parser = ml_sub.add_parser("predict", help="predict labels using a trained ViroFlow ML model")
    predict_parser.add_argument("--model", type=Path, required=True)
    predict_parser.add_argument("--features", type=Path, required=True)
    predict_parser.add_argument("--output", type=Path, required=True)
    predict_parser.add_argument("--sample-column", default="sample")

    anomaly_parser = ml_sub.add_parser("anomaly", help="score unusual genomes in a cohort")
    anomaly_parser.add_argument("--features", type=Path, required=True)
    anomaly_parser.add_argument("--output", type=Path, required=True)
    anomaly_parser.add_argument("--sample-column", default="sample")
    anomaly_parser.add_argument("--contamination", default="auto")
    anomaly_parser.add_argument(
        "--exclude-column",
        action="append",
        default=[],
        help="non-feature column to ignore; may be provided multiple times",
    )
    return parser


def _init_project(directory: Path, preset: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    files = {
        "config.yaml": SARS_COV_2_CONFIG_TEMPLATE if preset == "sars-cov-2" else CONFIG_TEMPLATE,
        "samples.csv": SAMPLES_TEMPLATE,
        "profile.yaml": PROFILE_TEMPLATE,
    }
    existing = [name for name in files if (directory / name).exists()]
    if existing:
        raise ConfigError("refusing to overwrite existing file(s): " + ", ".join(existing))
    for name, content in files.items():
        (directory / name).write_text(content, encoding="utf-8", newline="\n")
    for child in ("reads", "references", "datasets", "models"):
        (directory / child).mkdir(exist_ok=True)
    print(f"Initialized {preset} ViroFlow project in {directory.resolve()}")


def _parse_contamination(value: str) -> str | float:
    if value == "auto":
        return value
    parsed = float(value)
    if not 0 < parsed <= 0.5:
        raise MLError("--contamination must be 'auto' or a float in (0, 0.5]")
    return parsed


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    print(banner())
    try:
        if args.command == "init":
            _init_project(args.directory, args.preset)
        elif args.command == "validate":
            config = load_yaml(args.config.resolve())
            advisories = validate_config(config, check_files=not args.no_file_check)
            print("Configuration is valid.")
            for advisory in advisories:
                print(f"Advisory: {advisory}")
        elif args.command == "run":
            completed = run_pipeline(args.config.resolve(), dry_run=args.dry_run, resume=args.resume)
            if completed:
                summary = completed[0].parent / "run_summary.tsv"
                write_summary(completed, summary)
                print(f"Completed {len(completed)} sample(s). Summary: {summary}")
        elif args.command == "analyze":
            profile = load_yaml(args.profile.resolve())
            result = analyze_genome(args.sample, args.input.resolve(), profile)
            paths = write_reports(result, args.output.resolve())
            print("Reports:")
            for name, path in paths.items():
                print(f"  {name}: {path}")
        elif args.command == "doctor":
            if args.config:
                config = load_yaml(args.config.resolve())
                names = required_tools_for_config(config)
                status = {name: tool_status().get(name) for name in names}
            else:
                status = tool_status()
            if args.as_json:
                print(json.dumps(status, indent=2))
            else:
                for name, path in status.items():
                    print(f"{name:12} {'OK ' + path if path else 'MISSING'}")
            return 0 if all(status.values()) else 1
        elif args.command == "ml":
            if args.ml_command == "features":
                output = extract_features(args.results.resolve(), args.output.resolve())
                print(f"Feature table: {output}")
            elif args.ml_command == "train":
                metrics = train_classifier(
                    args.features.resolve(),
                    args.label_column,
                    args.model.resolve(),
                    args.report.resolve(),
                    sample_column=args.sample_column,
                )
                print("Model trained. Cross-validation metrics:")
                print(json.dumps(metrics, indent=2))
            elif args.ml_command == "predict":
                output = predict_classifier(
                    args.model.resolve(),
                    args.features.resolve(),
                    args.output.resolve(),
                    sample_column=args.sample_column,
                )
                print(f"Predictions: {output}")
            elif args.ml_command == "anomaly":
                output = score_anomalies(
                    args.features.resolve(),
                    args.output.resolve(),
                    contamination=_parse_contamination(args.contamination),
                    sample_column=args.sample_column,
                    exclude_columns=set(args.exclude_column),
                )
                print(f"Anomaly scores: {output}")
        return 0
    except (ConfigError, MLError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
