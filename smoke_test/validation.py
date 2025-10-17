"""Output and catalog validation helpers."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .constants import CATALOG_ABS_TOL, CATALOG_REL_TOL, REPO_ROOT
from .errors import SmokeTestError


def verify_outputs(output_dir: Path) -> List[Path]:
    out_files = sorted(output_dir.glob("*.out"))
    lc_files = list(output_dir.rglob("*.lc"))
    if not out_files:
        raise SmokeTestError(f"No .out files found in {output_dir}")
    if not lc_files:
        raise SmokeTestError(f"No .lc files found under {output_dir}")
    if out_files[0].stat().st_size == 0:
        raise SmokeTestError(f"Summary file is empty: {out_files[0]}")
    return out_files


def _resolve_param_path(raw: str, role: str) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = REPO_ROOT / raw
    resolved = path.resolve()
    if role.endswith("directory") and not resolved.is_dir():
        raise SmokeTestError(f"{role.capitalize()} '{resolved}' does not exist or is not a directory")
    return resolved


def _load_catalog_paths(directory: Path, list_file: str, role: str) -> List[Path]:
    list_path = (directory / list_file).resolve()
    lines = list_path.read_text(encoding="utf-8").splitlines()

    catalog_paths: List[Path] = []
    for raw in lines:
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        catalog_file = parts[-1]
        catalog_path = (directory / catalog_file).resolve()
        if not catalog_path.is_file():
            raise SmokeTestError(
                f"{role.capitalize()} catalog '{catalog_file}' referenced in {list_path} not found at {catalog_path}"
            )
        catalog_paths.append(catalog_path)

    if not catalog_paths:
        raise SmokeTestError(f"No {role} catalog entries found in {list_path}")
    return catalog_paths


def _load_mass_dist_pairs(path: Path, role: str) -> List[Tuple[float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise SmokeTestError(f"{role.capitalize()} catalog {path} is empty")

    header = lines[0].strip().split()
    if not header:
        raise SmokeTestError(f"{role.capitalize()} catalog {path} has an empty header")

    mass_idx = header.index("Mass")
    dist_idx = header.index("Dist")

    pairs: List[Tuple[float, float]] = []
    max_idx = max(mass_idx, dist_idx)
    for line_number, raw in enumerate(lines[1:], start=2):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) <= max_idx:
            raise SmokeTestError(
                f"{path}:{line_number} has insufficient columns to read Mass/Dist (expected index {max_idx})"
            )
        mass = float(parts[mass_idx])
        dist = float(parts[dist_idx])
        pairs.append((mass, dist))

    if not pairs:
        raise SmokeTestError(f"{role.capitalize()} catalog {path} does not contain any data rows")
    return pairs


def _gather_catalog_pairs(
    params: Dict[str, str],
    dir_key: str,
    list_key: str,
    role: str,
) -> Tuple[List[Path], List[Tuple[float, float]]]:
    directory_value = params.get(dir_key)
    list_value = params.get(list_key)
    if not directory_value or not list_value:
        raise SmokeTestError(f"Parameter file is missing {dir_key} or {list_key} required for {role} catalogs")

    directory = _resolve_param_path(directory_value, f"{role} directory")
    catalog_paths = _load_catalog_paths(directory, list_value, role)
    pairs: List[Tuple[float, float]] = []
    for catalog_path in catalog_paths:
        pairs.extend(_load_mass_dist_pairs(catalog_path, role))
    return catalog_paths, pairs


def _pair_matches(target: Tuple[float, float], catalog: List[Tuple[float, float]]) -> bool:
    target_mass, target_dist = target
    if math.isnan(target_mass) or math.isnan(target_dist):
        return False
    for mass, dist in catalog:
        if math.isclose(
            target_mass,
            mass,
            rel_tol=CATALOG_REL_TOL,
            abs_tol=CATALOG_ABS_TOL,
        ) and math.isclose(
            target_dist,
            dist,
            rel_tol=CATALOG_REL_TOL,
            abs_tol=CATALOG_ABS_TOL,
        ):
            return True
    return False


def verify_catalog_alignment(out_files: Sequence[Path], params: Dict[str, str]) -> None:
    source_catalog_paths, source_pairs = _gather_catalog_pairs(params, "SOURCE_DIR", "SOURCE_LIST", "source")
    lens_catalog_paths, lens_pairs = _gather_catalog_pairs(params, "LENS_DIR", "LENS_LIST", "lens")

    lens_catalog_names = ", ".join(path.name for path in lens_catalog_paths)
    source_catalog_names = ", ".join(path.name for path in source_catalog_paths)

    for out_file in out_files:
        with out_file.open(encoding="utf-8") as handle:
            header_line = next(handle)

            header = header_line.strip().split()
            if not header:
                raise SmokeTestError(f"Output file {out_file} has an empty header row")

            source_mass_idx = header.index("Source_Mass")
            source_dist_idx = header.index("Source_Dist")
            lens_mass_idx = header.index("Lens_Mass")
            lens_dist_idx = header.index("Lens_Dist")

            max_idx = max(source_mass_idx, source_dist_idx, lens_mass_idx, lens_dist_idx)
            for row_number, raw in enumerate(handle, start=2):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) <= max_idx:
                    raise SmokeTestError(
                        f"{out_file}:{row_number} has insufficient columns to read Source/Lens Mass or Dist"
                    )

                src_mass = float(parts[source_mass_idx])
                src_dist = float(parts[source_dist_idx])
                lens_mass = float(parts[lens_mass_idx])
                lens_dist = float(parts[lens_dist_idx])

                if not _pair_matches((src_mass, src_dist), source_pairs):
                    raise SmokeTestError(
                        f"{out_file}:{row_number} Source mass/dist ({src_mass:.6g}, {src_dist:.6g}) "
                        f"not found in source catalogs: {source_catalog_names}"
                    )
                if not _pair_matches((lens_mass, lens_dist), lens_pairs):
                    raise SmokeTestError(
                        f"{out_file}:{row_number} Lens mass/dist ({lens_mass:.6g}, {lens_dist:.6g}) "
                        f"not found in lens catalogs: {lens_catalog_names}"
                    )


def verify_binary_source_columns(params: Dict[str, str]) -> None:
    """Verify source catalogs contain required columns when MULTIPLE_SOURCES=1."""
    multiple_sources = params.get("MULTIPLE_SOURCES", "0").strip()
    
    # Only validate if multiple sources is enabled
    if multiple_sources not in ("1", "1.0"):
        return
    
    # Required columns for binary source simulations
    # These are accessed via datadict in buildEvent.cpp
    required_named_cols = {"Is_Binary", "ID", "primary_ID", "combined_logP"}
    
    # These columns must exist by position (standard catalog columns)
    # They're accessed via Sources->RADIUS, Sources->DIST, etc.
    required_standard_cols = {"mul", "mub", "Mass", "Radius", "Dist"}
    
    directory_value = params.get("SOURCE_DIR")
    list_value = params.get("SOURCE_LIST")
    if not directory_value or not list_value:
        raise SmokeTestError(
            "MULTIPLE_SOURCES=1 requires SOURCE_DIR and SOURCE_LIST to be set"
        )
    
    directory = _resolve_param_path(directory_value, "source directory")
    catalog_paths = _load_catalog_paths(directory, list_value, "source")
    
    # Check each source catalog
    for catalog_path in catalog_paths:
        lines = catalog_path.read_text(encoding="utf-8").splitlines()
        
        if not lines:
            raise SmokeTestError(
                f"MULTIPLE_SOURCES=1 but source catalog {catalog_path} is empty"
            )
        
        # Find the header line (first non-comment, non-empty line)
        header_line = None
        for line in lines:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                header_line = stripped
                break
        
        if not header_line:
            raise SmokeTestError(
                f"MULTIPLE_SOURCES=1 but source catalog {catalog_path} has no header"
            )
        
        header_cols = set(header_line.split())
        
        # Check for required named columns
        missing_named = required_named_cols - header_cols
        if missing_named:
            raise SmokeTestError(
                f"MULTIPLE_SOURCES=1 but source catalog {catalog_path.name} "
                f"is missing required binary source columns: {', '.join(sorted(missing_named))}"
            )
        
        # Check for required standard columns
        missing_standard = required_standard_cols - header_cols
        if missing_standard:
            raise SmokeTestError(
                f"MULTIPLE_SOURCES=1 but source catalog {catalog_path.name} "
                f"is missing required standard columns: {', '.join(sorted(missing_standard))}"
            )


__all__ = ["verify_binary_source_columns", "verify_catalog_alignment", "verify_outputs"]
