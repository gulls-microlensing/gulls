"""Helpers for invoking gulls binaries."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, Sequence


def run_command(
    cmd: Sequence[str],
    env: Dict[str, str],
    timeout: float | None,
    *,
    cwd: Path,
) -> subprocess.CompletedProcess[str]:
    effective_timeout = None if timeout is None or timeout <= 0 else timeout
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=effective_timeout,
    )


__all__ = ["run_command"]
