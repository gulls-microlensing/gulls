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


def verify_source_lens_compatibility(params: Dict[str, str]) -> None:
    """Verify that at least some valid source/lens pairs exist in catalogs."""
    _, source_pairs = _gather_catalog_pairs(params, "SOURCE_DIR", "SOURCE_LIST", "source")
    _, lens_pairs = _gather_catalog_pairs(params, "LENS_DIR", "LENS_LIST", "lens")
    
    # Check 1: At least one valid pair where source_dist > lens_dist
    valid_pairs = 0
    for src_mass, src_dist in source_pairs:
        for lens_mass, lens_dist in lens_pairs:
            if src_dist > lens_dist:
                valid_pairs += 1
                if valid_pairs >= 10:  # Early exit after finding enough
                    break
        if valid_pairs >= 10:
            break
    
    if valid_pairs == 0:
        src_dists = [d for _, d in source_pairs]
        lens_dists = [d for _, d in lens_pairs]
        raise SmokeTestError(
            f"No valid source/lens pairs found!\n"
            f"  All sources must be farther than at least some lenses (source_dist > lens_dist).\n"
            f"  Source distances: min={min(src_dists):.2f}, max={max(src_dists):.2f}, median={sorted(src_dists)[len(src_dists)//2]:.2f} kpc\n"
            f"  Lens distances: min={min(lens_dists):.2f}, max={max(lens_dists):.2f}, median={sorted(lens_dists)[len(lens_dists)//2]:.2f} kpc\n"
            f"  Check catalog generation - sources should generally be farther than lenses."
        )
    
    # Check 2: Verify distances are reasonable (not NaN/Inf, positive, within galaxy scale)
    for role, pairs in [("source", source_pairs), ("lens", lens_pairs)]:
        for mass, dist in pairs:
            if math.isnan(dist) or math.isinf(dist):
                raise SmokeTestError(f"{role.capitalize()} catalog contains NaN or Inf distance values")
            if dist <= 0:
                raise SmokeTestError(f"{role.capitalize()} catalog contains non-positive distance: {dist} kpc")
            if dist > 50:  # Milky Way is ~50 kpc diameter
                raise SmokeTestError(f"{role.capitalize()} catalog contains implausibly large distance: {dist} kpc (>50 kpc)")


def verify_catalog_columns(params: Dict[str, str]) -> None:
    """Verify source and lens catalogs contain all required standard columns."""
    # Required columns for all source catalogs
    required_source_cols = {"mul", "mub", "Mass", "Radius", "Dist"}
    
    # Required columns for all lens catalogs
    required_lens_cols = {"mul", "mub", "Mass", "Dist"}
    
    # Validate source catalogs
    directory_value = params.get("SOURCE_DIR")
    list_value = params.get("SOURCE_LIST")
    if directory_value and list_value:
        directory = _resolve_param_path(directory_value, "source directory")
        catalog_paths = _load_catalog_paths(directory, list_value, "source")
        
        for catalog_path in catalog_paths:
            lines = catalog_path.read_text(encoding="utf-8").splitlines()
            
            if not lines:
                raise SmokeTestError(f"Source catalog {catalog_path} is empty")
            
            # Find the header line (first non-comment, non-empty line)
            header_line = None
            for line in lines:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    header_line = stripped
                    break
            
            if not header_line:
                raise SmokeTestError(f"Source catalog {catalog_path} has no header")
            
            header_cols = set(header_line.split())
            missing = required_source_cols - header_cols
            if missing:
                raise SmokeTestError(
                    f"Source catalog {catalog_path.name} is missing required columns: {', '.join(sorted(missing))}"
                )
    
    # Validate lens catalogs
    directory_value = params.get("LENS_DIR")
    list_value = params.get("LENS_LIST")
    if directory_value and list_value:
        directory = _resolve_param_path(directory_value, "lens directory")
        catalog_paths = _load_catalog_paths(directory, list_value, "lens")
        
        for catalog_path in catalog_paths:
            lines = catalog_path.read_text(encoding="utf-8").splitlines()
            
            if not lines:
                raise SmokeTestError(f"Lens catalog {catalog_path} is empty")
            
            # Find the header line (first non-comment, non-empty line)
            header_line = None
            for line in lines:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    header_line = stripped
                    break
            
            if not header_line:
                raise SmokeTestError(f"Lens catalog {catalog_path} has no header")
            
            header_cols = set(header_line.split())
            missing = required_lens_cols - header_cols
            if missing:
                raise SmokeTestError(
                    f"Lens catalog {catalog_path.name} is missing required columns: {', '.join(sorted(missing))}"
                )


def verify_binary_source_columns(params: Dict[str, str]) -> None:
    """Verify source catalogs contain required columns when MULTIPLE_SOURCES=1."""
    multiple_sources = params.get("MULTIPLE_SOURCES", "0").strip()
    
    # Only validate if multiple sources is enabled
    if multiple_sources not in ("1", "1.0"):
        return
    
    # Required columns for binary source simulations (beyond standard columns)
    # These are accessed via datadict in buildEvent.cpp
    required_binary_cols = {"Is_Binary", "ID", "primary_ID", "combined_logP"}
    
    directory_value = params.get("SOURCE_DIR")
    list_value = params.get("SOURCE_LIST")
    if not directory_value or not list_value:
        raise SmokeTestError(
            "MULTIPLE_SOURCES=1 requires SOURCE_DIR and SOURCE_LIST to be set"
        )
    
    directory = _resolve_param_path(directory_value, "source directory")
    catalog_paths = _load_catalog_paths(directory, list_value, "source")
    
    # Check each source catalog for binary-specific columns
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
        
        # Check for required binary-specific columns
        missing = required_binary_cols - header_cols
        if missing:
            raise SmokeTestError(
                f"MULTIPLE_SOURCES=1 but source catalog {catalog_path.name} "
                f"is missing required binary source columns: {', '.join(sorted(missing))}"
            )


def verify_weather_coverage(params: Dict[str, str]) -> None:
    """Verify weather file covers the full simulation duration."""
    sim_length_str = params.get("SIMULATION_LENGTH")
    weather_dir = params.get("WEATHER_PROFILE_DIR")
    
    if not sim_length_str or not weather_dir:
        # If not specified, we can't validate
        return
    
    try:
        sim_length = float(sim_length_str)
    except ValueError:
        raise SmokeTestError(f"SIMULATION_LENGTH must be a number, got: {sim_length_str}")
    
    # Weather files are referenced in observatory files, but we can check the base directory
    weather_dir_path = _resolve_param_path(weather_dir, "weather directory")
    
    # Look for .weather files in the directory
    weather_files = list(weather_dir_path.glob("*.weather"))
    
    if not weather_files:
        # No weather files found, might be okay if observatories don't use them
        return
    
    # Check the first weather file we find
    for weather_file in weather_files:
        lines = weather_file.read_text(encoding="utf-8").splitlines()
        
        if not lines:
            raise SmokeTestError(f"Weather file {weather_file.name} is empty")
        
        # Parse weather data (format: time observing_flag)
        max_time = 0.0
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    time = float(parts[0])
                    max_time = max(max_time, time)
                except ValueError:
                    continue
        
        if max_time < sim_length:
            raise SmokeTestError(
                f"Weather file {weather_file.name} only covers {max_time:.2f} days "
                f"but SIMULATION_LENGTH is {sim_length:.2f} days.\n"
                f"  Weather files must cover the entire simulation duration.\n"
                f"  Use scripts/make_weather.py to generate appropriate weather profiles."
            )
        
        # Only check the first file
        break


def verify_rates_file(params: Dict[str, str]) -> None:
    """Verify rates file has valid parameter ranges."""
    rates_file_str = params.get("RATES_FILE")
    
    if not rates_file_str:
        # Rates file is optional for some simulation types
        return
    
    rates_file = _resolve_param_path(rates_file_str, "rates file")
    
    if not rates_file.exists():
        raise SmokeTestError(f"Rates file not found: {rates_file}")
    
    lines = rates_file.read_text(encoding="utf-8").splitlines()
    
    if not lines:
        raise SmokeTestError(f"Rates file {rates_file.name} is empty")
    
    # Parse rates data (format: field u0min u0max t0min t0max tEmin tEmax rate)
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        
        parts = line.split()
        if len(parts) < 8:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: expected 8 columns "
                f"(field u0min u0max t0min t0max tEmin tEmax rate), got {len(parts)}"
            )
        
        try:
            field, u0min, u0max, t0min, t0max, tEmin, tEmax, rate = [float(x) for x in parts[:8]]
        except ValueError as e:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: invalid numeric value: {e}"
            )
        
        # Validate ranges
        if u0min >= u0max:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: u0min ({u0min}) >= u0max ({u0max})"
            )
        
        if t0min >= t0max:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: t0min ({t0min}) >= t0max ({t0max})"
            )
        
        if tEmin >= tEmax:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: tEmin ({tEmin}) >= tEmax ({tEmax})"
            )
        
        if rate < 0:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: rate ({rate}) is negative"
            )
        
        if u0min < 0:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: u0min ({u0min}) is negative"
            )
        
        if tEmin <= 0:
            raise SmokeTestError(
                f"Rates file {rates_file.name} line {line_num}: tEmin ({tEmin}) must be positive"
            )


def verify_nfilters_matches_catalogs(params: Dict[str, str]) -> None:
    """Verify NFILTERS parameter matches the number of magnitude columns in catalogs."""
    nfilters_str = params.get("NFILTERS")
    
    if not nfilters_str:
        # NFILTERS not specified, can't validate
        return
    
    try:
        nfilters = int(nfilters_str)
    except ValueError:
        raise SmokeTestError(f"NFILTERS must be an integer, got: {nfilters_str}")
    
    # Check source catalogs
    directory_value = params.get("SOURCE_DIR")
    list_value = params.get("SOURCE_LIST")
    
    if directory_value and list_value:
        directory = _resolve_param_path(directory_value, "source directory")
        catalog_paths = _load_catalog_paths(directory, list_value, "source")
        
        for catalog_path in catalog_paths[:1]:  # Just check the first one
            lines = catalog_path.read_text(encoding="utf-8").splitlines()
            
            if not lines:
                continue
            
            # Find header line
            header_line = None
            for line in lines:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    header_line = stripped
                    break
            
            if not header_line:
                continue
            
            header_cols = header_line.split()
            total_cols = len(header_cols)
            
            # The first NFILTERS columns should be magnitude columns
            # Common non-magnitude columns that should come after magnitudes
            physical_cols = {"mul", "mub", "Mass", "Radius", "Dist", "Vr", "U", "V", "W", 
                           "RA2000.0", "DEC2000.0", "Teff", "logg", "[Fe/H]", "l", "b",
                           "x", "y", "z", "pop", "age", "CL", "Av", "[alpha/Fe]", "Mbol",
                           "iMass", "In_Final_Phase", "Fe/H_initial", "Fe/H_evolved",
                           "Is_Binary", "ID", "primary_ID", "combined_logP", "combined_logL",
                           "total_mass", "q", "logL", "logTeff", "log_radius", "vr_bc", "VR_LSR"}
            
            # Count columns before we hit physical property columns
            magnitude_cols = 0
            for col in header_cols:
                if col in physical_cols:
                    break
                magnitude_cols += 1
            
            if magnitude_cols != nfilters:
                raise SmokeTestError(
                    f"NFILTERS={nfilters} but source catalog {catalog_path.name} has "
                    f"{magnitude_cols} magnitude columns before physical properties.\n"
                    f"  The first {nfilters} columns must be photometric magnitudes.\n"
                    f"  Check that NFILTERS matches your catalog format."
                )
            
            # Only check the first catalog
            break


__all__ = [
    "verify_binary_source_columns",
    "verify_catalog_alignment", 
    "verify_catalog_columns",
    "verify_nfilters_matches_catalogs",
    "verify_outputs",
    "verify_rates_file",
    "verify_source_lens_compatibility",
    "verify_weather_coverage",
]
