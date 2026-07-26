from __future__ import annotations

import json
import shlex
import subprocess
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path


def display_command(command: Sequence[str]) -> str:
    return shlex.join(str(part) for part in command)


class Runner:
    def __init__(self, log_dir: Path, dry_run: bool = False):
        self.log_dir = log_dir
        self.dry_run = dry_run
        self.commands: list[dict[str, object]] = []
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        command: Sequence[str | Path],
        label: str,
        stdout_path: Path | None = None,
    ) -> None:
        normalized = [str(part) for part in command]
        self._record(label, [normalized])
        if self.dry_run:
            print(f"[dry-run] {display_command(normalized)}")
            return
        log_path = self.log_dir / f"{label}.log"
        with log_path.open("w", encoding="utf-8") as log_handle:
            if stdout_path:
                stdout_path.parent.mkdir(parents=True, exist_ok=True)
                with stdout_path.open("wb") as output_handle:
                    subprocess.run(
                        normalized,
                        check=True,
                        stdout=output_handle,
                        stderr=log_handle,
                    )
            else:
                subprocess.run(
                    normalized,
                    check=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                )

    def pipe(self, commands: Sequence[Sequence[str | Path]], label: str) -> None:
        normalized = [[str(part) for part in command] for command in commands]
        self._record(label, normalized)
        if self.dry_run:
            print("[dry-run] " + " | ".join(display_command(command) for command in normalized))
            return
        log_path = self.log_dir / f"{label}.log"
        processes: list[subprocess.Popen[bytes]] = []
        with log_path.open("wb") as log_handle:
            previous_stdout = None
            for index, command in enumerate(normalized):
                process = subprocess.Popen(
                    command,
                    stdin=previous_stdout,
                    stdout=subprocess.PIPE if index < len(normalized) - 1 else log_handle,
                    stderr=log_handle,
                )
                if previous_stdout is not None:
                    previous_stdout.close()
                previous_stdout = process.stdout
                processes.append(process)
            return_codes = [process.wait() for process in processes]
        failures = [
            f"{display_command(command)} (exit {code})"
            for command, code in zip(normalized, return_codes)
            if code != 0
        ]
        if failures:
            raise subprocess.CalledProcessError(
                return_codes[-1],
                normalized[-1],
                output="; ".join(failures),
            )

    def write_manifest(self, path: Path, metadata: dict[str, object]) -> None:
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "dry_run": self.dry_run,
            "metadata": metadata,
            "commands": self.commands,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    def _record(self, label: str, commands: list[list[str]]) -> None:
        self.commands.append({"label": label, "commands": commands})

