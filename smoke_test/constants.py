"""Static configuration shared across smoke test modules."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple


REPO_ROOT: Path = Path(__file__).resolve().parents[1]
BUILD_BIN_DEFAULT: Path = REPO_ROOT / "bin"
PARAM_DIR: Path = REPO_ROOT / "smoke_test" / "parameterfiles"
CATALOG_REL_TOL: float = 1e-9
CATALOG_ABS_TOL: float = 1e-8

CaseDef = Tuple[str, str, str]

CASES: Tuple[CaseDef, ...] = (
    ("std-single", "gulls_std.x", "smoke_std.prm"),
    ("std-binary", "gulls_std.x", "smoke_std_binary.prm"),
    ("std-heavy", "gulls_std.x", "smoke_std_heavy.prm"),
    ("croin-single", "gulls_croin.x", "smoke_croin.prm"),
    ("croin-binary", "gulls_croin.x", "smoke_croin_binary.prm"),
    ("croin-heavy", "gulls_croin.x", "smoke_croin_heavy.prm"),
    ("fish-single", "gullsFish.x", "smoke_fish.prm"),
    ("fish-binary", "gullsFish.x", "smoke_fish_binary.prm"),
    ("fish-heavy", "gullsFish.x", "smoke_fish_heavy.prm"),
    # Houston catalog tests - different seeds to test for serendipitous success
    ("std-houston-seed1", "gulls_std.x", "smoke_std_houston_seed1.prm"),
    ("std-houston-seed2", "gulls_std.x", "smoke_std_houston_seed2.prm"),
    ("std-houston-seed3", "gulls_std.x", "smoke_std_houston_seed3.prm"),
)

CASE_LABELS = {label for label, _, _ in CASES}
CASE_LOOKUP: Dict[str, CaseDef] = {label: (label, exe, prm) for label, exe, prm in CASES}
CASE_EXEC_MAP: Dict[str, List[CaseDef]] = {}
for case in CASES:
    CASE_EXEC_MAP.setdefault(case[1], []).append(case)

CASE_CHOICES: Tuple[str, ...] = tuple(sorted(CASE_LABELS | set(CASE_EXEC_MAP)))

__all__ = [
    "BUILD_BIN_DEFAULT",
    "CASE_CHOICES",
    "CASE_EXEC_MAP",
    "CASE_LABELS",
    "CASE_LOOKUP",
    "CASES",
    "CaseDef",
    "CATALOG_ABS_TOL",
    "CATALOG_REL_TOL",
    "PARAM_DIR",
    "REPO_ROOT",
]
