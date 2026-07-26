from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .analysis import analyze_genome
from .config import ConfigError, load_yaml
from .pipeline import run_pipeline, tool_status, write_summary
from .report import write_reports
from .templates import CONFIG_TEMPLATE, PROFILE_TEMPLATE, SAMPLES_TEMPLATE


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viroflow",
        description="Viral raw-read workflow and comparative antigenic evidence screen.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="write starter configuration files")
    init_parser.add_argument("directory", nargs="?", type=Path, default=Path("viroflow-project"))

    run_parser = subparsers.add_parser("run", help="run raw-read workflow")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--dry-run", action="store_true", help="print commands without executing")

    analyze_parser = subparsers.add_parser(
        "analyze", help="analyze an assembled genome FASTA with an evidence profile"
    )
    analyze_parser.add_argument("--profile", type=Path, required=True)
    analyze_parser.add_argument("--input", type=Path, required=True)
    analyze_parser.add_argument("--sample", default="sample")
    analyze_parser.add_argument("--output", type=Path, required=True)

    doctor_parser = subparsers.add_parser("doctor", help="check external bioinformatics tools")
    doctor_parser.add_argument("--json", action="store_true", dest="as_json")
    return parser


def _init_project(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    files = {
        "config.yaml": CONFIG_TEMPLATE,
        "samples.csv": SAMPLES_TEMPLATE,
        "profile.yaml": PROFILE_TEMPLATE,
    }
    existing = [name for name in files if (directory / name).exists()]
    if existing:
        raise ConfigError("refusing to overwrite existing file(s): " + ", ".join(existing))
    for name, content in files.items():
        (directory / name).write_text(content, encoding="utf-8", newline="\n")
    for child in ("reads", "references", "datasets"):
        (directory / child).mkdir(exist_ok=True)
    print(f"Initialized ViroFlow project in {directory.resolve()}")


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "init":
            _init_project(args.directory)
        elif args.command == "run":
            completed = run_pipeline(args.config.resolve(), dry_run=args.dry_run)
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
            status = tool_status()
            if args.as_json:
                print(json.dumps(status, indent=2))
            else:
                for name, path in status.items():
                    print(f"{name:12} {'OK ' + path if path else 'MISSING'}")
            return 0 if all(status.values()) else 1
        return 0
    except (ConfigError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

