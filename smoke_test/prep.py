"""Case preparation utilities for the smoke test."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .constants import PARAM_DIR, REPO_ROOT, CaseDef
from .errors import SmokeTestError


@dataclass
class PreparedCase:
    """Executable invocation together with resolved parameter metadata."""

    label: str
    exe_name: str
    exe_path: Path
    param_path: Path
    exec_param_path: Path
    run_name: str
    output_root: Path
    output_dir: Path
    params: Dict[str, str]


def ensure_executable(path: Path) -> None:
    if not path.is_file():
        raise SmokeTestError(f"Missing executable: {path}")
    if not os.access(path, os.X_OK):
        raise SmokeTestError(f"Executable is not runnable: {path}")


def parse_parameter_file(path: Path) -> Dict[str, str]:
    params: Dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip()
    return params


def _discover_weather_file(params: Dict[str, str]) -> tuple[str | None, str | None]:
    weather_dir = params.get("WEATHER_PROFILE_DIR")
    weather_file = params.get("WEATHER_PROFILE")
    if weather_dir and weather_file:
        return weather_dir, weather_file

    obs_dir = params.get("OBSERVATORY_DIR")
    obs_list = params.get("OBSERVATORY_LIST")
    if not (weather_dir and obs_dir and obs_list):
        return weather_dir, weather_file

    list_path = (REPO_ROOT / obs_dir / obs_list).resolve()
    if not list_path.is_file():
        return weather_dir, weather_file

    entries = [
        line.strip()
        for line in list_path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]

    for entry in entries:
        obs_path = (REPO_ROOT / obs_dir / entry).resolve()
        if not obs_path.is_file():
            continue
        for raw in obs_path.read_text(encoding="utf-8").splitlines():
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) >= 2 and parts[0].upper() == "WEATHER_PROFILE":
                return weather_dir, parts[1]

    return weather_dir, weather_file


def ensure_weather_profile(params: Dict[str, str]) -> None:
    weather_dir, weather_file = _discover_weather_file(params)
    num_days_raw = params.get("NUM_SIM_DAYS")

    if not weather_dir or not weather_file or not num_days_raw:
        return

    num_days = int(float(num_days_raw))

    if num_days <= 0:
        return

    weather_path = (REPO_ROOT / weather_dir / weather_file).resolve()
    required_entries = int(round((num_days + 1) * 4))

    lines: List[str] = weather_path.read_text(encoding="utf-8").splitlines() if weather_path.is_file() else []

    values: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) >= 2:
            values.append(float(parts[1]))

    if not values:
        values = [1.0]

    pattern = list(values)
    expanded = [pattern[idx % len(pattern)] for idx in range(required_entries)]

    weather_path.parent.mkdir(parents=True, exist_ok=True)
    with weather_path.open("w", encoding="utf-8") as handle:
        for idx, value in enumerate(expanded):
            handle.write(f"{idx / 4:.2f} {value:.6g}\n")


def _write_param_override(src: Path, dest: Path, run_name: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    replaced = False
    with src.open(encoding="utf-8") as reader, dest.open("w", encoding="utf-8") as writer:
        for line in reader:
            if line.strip().startswith("RUN_NAME="):
                writer.write(f"RUN_NAME={run_name}\n")
                replaced = True
            else:
                writer.write(line)
        if not replaced:
            writer.write(f"\nRUN_NAME={run_name}\n")


def prepare_cases(build_bin: Path, selected: Sequence[CaseDef]) -> Tuple[List[PreparedCase], List[str]]:
    prepared: List[PreparedCase] = []
    failures: List[str] = []
    run_name_usage: Dict[str, int] = {}

    for label, exe_name, prm_filename in selected:
        exe_path = build_bin / exe_name
        if not exe_path.is_file():
            failures.append(f"{label}: Missing executable: {exe_path}")
            continue
        if not os.access(exe_path, os.X_OK):
            failures.append(f"{label}: Executable is not runnable: {exe_path}")
            continue

        param_path = PARAM_DIR / prm_filename
        if not param_path.is_file():
            failures.append(f"{label}: Parameter file not found: {param_path}")
            continue

        params = parse_parameter_file(param_path)
        ensure_weather_profile(params)
        base_run_name = params.get("RUN_NAME", param_path.stem)
        usage_count = run_name_usage.get(base_run_name, 0)
        run_name_usage[base_run_name] = usage_count + 1
        if usage_count == 0:
            run_name = base_run_name
            exec_param_path = param_path
        else:
            run_name = f"{base_run_name}_{label}"
            params = params.copy()
            params["RUN_NAME"] = run_name
            exec_param_path = PARAM_DIR / "__generated__" / f"{param_path.stem}_{label}.prm"
            _write_param_override(param_path, exec_param_path, run_name)

        output_dir_value = params.get("OUTPUT_DIR")
        if not output_dir_value:
            failures.append(f"{label}: OUTPUT_DIR missing in {param_path}")
            continue

        output_root = (REPO_ROOT / output_dir_value).resolve()
        output_dir = output_root / run_name
        prepared.append(
            PreparedCase(
                label=label,
                exe_name=exe_name,
                exe_path=exe_path,
                param_path=param_path,
                exec_param_path=exec_param_path,
                run_name=run_name,
                output_root=output_root,
                output_dir=output_dir,
                params=params,
            )
        )

    return prepared, failures


__all__ = [
    "PreparedCase",
    "ensure_executable",
    "ensure_weather_profile",
    "parse_parameter_file",
    "prepare_cases",
]
