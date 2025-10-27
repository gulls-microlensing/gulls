"""Summary extraction utilities for smoke test outputs."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Sequence, Tuple

from .errors import SmokeTestError

SUMMARY_COLUMN_MAP = {
    "lens_mass": "Lens_Mass",
    "lens_dist": "Lens_Dist",
    "source_dist": "Source_Dist",
    "theta_e": "thetaE",
    "source_mul": "Source_mul",
    "source_mub": "Source_mub",
    "source_l": "galactic_l",
    "source_b": "galactic_b",
    "lens_mul": "Lens_mul",
    "lens_mub": "Lens_mub",
    "lens_l": "Lens_l",
    "lens_b": "Lens_b",
    "pi_n": "piEN",
    "pi_e": "piEE",
    "u0": "u0lens1",
    "alpha_event": "alpha",
    "t0": "t0lens1",
    "tE_ref": "tE_ref",
    "rho": "rho",
    "event_ra": "ra_deg",
    "event_dec": "dec_deg",
    "pm_helio_alpha": "murel_helio_alpha",
    "pm_helio_delta": "murel_helio_delta",
    "pm_ref_alpha": "murel_ref_alpha",
    "pm_ref_delta": "murel_ref_delta",
}

EVENT_ID_COLUMNS = ("EventID", "SubRun", "Field")


def _extract_summary_metrics(out_file: Path) -> Dict[Tuple[int, int, int], Dict[str, float]]:
    with out_file.open(encoding="utf-8") as handle:
        header_line = next(handle)

        header = header_line.strip().split()
        if not header:
            raise SmokeTestError(f"Output file {out_file} has an empty header row")

        event_idx = header.index(EVENT_ID_COLUMNS[0])
        subrun_idx = header.index(EVENT_ID_COLUMNS[1])
        field_idx = header.index(EVENT_ID_COLUMNS[2])

        col_indices: Dict[str, int | None] = {}
        for key, column_name in SUMMARY_COLUMN_MAP.items():
            col_indices[key] = header.index(column_name) if column_name in header else None

        required_indices = [idx for idx in col_indices.values() if idx is not None] + [
            event_idx,
            subrun_idx,
            field_idx,
        ]
        max_index = max(required_indices)

        metrics: Dict[Tuple[int, int, int], Dict[str, float]] = {}
        for row_number, raw in enumerate(handle, start=2):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) <= max_index:
                raise SmokeTestError(
                    f"{out_file}:{row_number} has insufficient columns to read summary metrics"
                )

            event_key = (
                int(float(parts[event_idx])),
                int(float(parts[subrun_idx])),
                int(float(parts[field_idx])),
            )

            summary: Dict[str, float] = {}
            for metric_key, col_idx in col_indices.items():
                if col_idx is None:
                    summary[metric_key] = math.nan
                else:
                    value_str = parts[col_idx]
                    summary[metric_key] = float(value_str)

            metrics[event_key] = summary

        if not metrics:
            raise SmokeTestError(f"No data rows found in {out_file}")

        return metrics


def gather_case_metrics(out_files: Sequence[Path]) -> Dict[Tuple[int, int, int], Dict[str, float]]:
    aggregated: Dict[Tuple[int, int, int], Dict[str, float]] = {}
    for out_file in out_files:
        aggregated.update(_extract_summary_metrics(out_file))
    return aggregated


__all__ = ["EVENT_ID_COLUMNS", "SUMMARY_COLUMN_MAP", "gather_case_metrics"]
