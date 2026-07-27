from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def collect_tool_versions(tools: list[str]) -> dict[str, dict[str, str | None]]:
    versions: dict[str, dict[str, str | None]] = {}
    for tool in tools:
        executable = shutil.which(tool)
        if not executable:
            versions[tool] = {"path": None, "version": None}
            continue
        try:
            completed = subprocess.run(
                [executable, "--version"],
                check=False,
                capture_output=True,
                text=True,
                timeout=15,
            )
            output = (completed.stdout or completed.stderr).strip().splitlines()
            version = output[0] if output else "version output unavailable"
        except (OSError, subprocess.TimeoutExpired) as exc:
            version = f"version check failed: {exc}"
        versions[tool] = {"path": executable, "version": version}
    return versions

