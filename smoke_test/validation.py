"""Output and catalog validation helpers."""
from __future__ import annotations

import hashlib
import math
import struct
import warnings
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .constants import CATALOG_ABS_TOL, CATALOG_REL_TOL, REPO_ROOT
from .errors import SmokeTestError


_CANONICAL_PSF_HASHES: Dict[str, str] = {
    # Canonical Roman-like PSF used by smoke tests (Nkern=145, Nsub=9, pixscale=0.11)
    "smoke_test/assets/observatories/WFI_PSF.psf": "ea9c8b2a2d3b0687487f8283661fd491fad492bb5e012776b8cf17324d65222b",
}


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


def _parse_key_value_file(path: Path) -> Dict[str, str]:
    """Parse simple key/value configuration files used for observatory and detector data."""
    result: Dict[str, str] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        key = parts[0]
        value = parts[1].strip() if len(parts) > 1 else ""
        result[key] = value
    return result


def _safe_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _safe_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _read_psf_header(psf_bytes: bytes) -> Tuple[int, Tuple[int, ...], Tuple[float, ...], int, int]:
    """Read the PSF binary header and return metadata."""
    offset = 0

    def take(fmt: str):
        nonlocal offset
        size = struct.calcsize(fmt)
        if offset + size > len(psf_bytes):
            raise SmokeTestError("PSF file truncated before header could be read")
        values = struct.unpack_from(fmt, psf_bytes, offset)
        offset += size
        return values

    version, = take("<i")
    nints, = take("<i")
    ivars = take(f"<{nints}i")
    ndoub, = take("<i")
    dvars = take(f"<{ndoub}d")
    nlong, = take("<i")
    if nlong:
        take(f"<{nlong}q")  # discard optional longs
    header_bytes = offset
    return version, ivars, dvars, header_bytes, len(psf_bytes)


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


def verify_input_files_exist(params: Dict[str, str]) -> None:
    """Verify all input files referenced in parameter file actually exist."""
    missing_files = []
    
    # Check observatory files
    obs_dir_str = params.get("OBSERVATORY_DIR")
    obs_list_str = params.get("OBSERVATORY_LIST")
    
    if obs_dir_str and obs_list_str:
        obs_dir = _resolve_param_path(obs_dir_str, "observatory directory")
        obs_list = _resolve_param_path(f"{obs_dir_str}/{obs_list_str}", "observatory list")
        
        if not obs_list.exists():
            missing_files.append(f"Observatory list: {obs_list}")
        else:
            # Check each observatory file listed
            lines = obs_list.read_text(encoding="utf-8").splitlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                obs_file = obs_dir / line
                if not obs_file.exists():
                    missing_files.append(f"Observatory file: {obs_file}")
                else:
                    # Check sequence file referenced in observatory
                    obs_content = obs_file.read_text(encoding="utf-8")
                    for obs_line in obs_content.splitlines():
                        obs_line = obs_line.strip()
                        if obs_line.startswith("OBSERVATION_SEQUENCE"):
                            parts = obs_line.split()
                            if len(parts) >= 2:
                                seq_file = obs_dir / parts[1]
                                if not seq_file.exists():
                                    missing_files.append(f"Sequence file: {seq_file} (referenced in {obs_file.name})")
    
    # Check source files
    source_dir_str = params.get("SOURCE_DIR")
    source_list_str = params.get("SOURCE_LIST")
    
    if source_dir_str and source_list_str:
        source_dir = _resolve_param_path(source_dir_str, "source directory")
        source_list = _resolve_param_path(f"{source_dir_str}/{source_list_str}", "source list")
        
        if not source_list.exists():
            missing_files.append(f"Source list: {source_list}")
        else:
            # Check each source catalog listed
            lines = source_list.read_text(encoding="utf-8").splitlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split()
                if len(parts) >= 6:  # field_number l b l_width b_width catalog_filename
                    catalog_file = source_dir / parts[5]
                    if not catalog_file.exists():
                        missing_files.append(f"Source catalog: {catalog_file}")
    
    # Check lens files
    lens_dir_str = params.get("LENS_DIR")
    lens_list_str = params.get("LENS_LIST")
    
    if lens_dir_str and lens_list_str:
        lens_dir = _resolve_param_path(lens_dir_str, "lens directory")
        lens_list = _resolve_param_path(f"{lens_dir_str}/{lens_list_str}", "lens list")
        
        if not lens_list.exists():
            missing_files.append(f"Lens list: {lens_list}")
        else:
            # Check each lens catalog listed
            lines = lens_list.read_text(encoding="utf-8").splitlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split()
                if len(parts) >= 6:  # field_number l b l_width b_width catalog_filename
                    catalog_file = lens_dir / parts[5]
                    if not catalog_file.exists():
                        missing_files.append(f"Lens catalog: {catalog_file}")
    
    # Check starfield files
    starfield_dir_str = params.get("STARFIELD_DIR")
    starfield_list_str = params.get("STARFIELD_LIST")
    
    if starfield_dir_str and starfield_list_str:
        starfield_dir = _resolve_param_path(starfield_dir_str, "starfield directory")
        starfield_list = _resolve_param_path(f"{starfield_dir_str}/{starfield_list_str}", "starfield list")
        
        if not starfield_list.exists():
            missing_files.append(f"Starfield list: {starfield_list}")
        else:
            # Check each starfield catalog listed
            lines = starfield_list.read_text(encoding="utf-8").splitlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split()
                if len(parts) >= 4:  # field_number l b area_or_weight catalog_filename
                    catalog_file = starfield_dir / parts[3]
                    if not catalog_file.exists():
                        missing_files.append(f"Starfield catalog: {catalog_file}")
    
    # Check weather files
    weather_dir_str = params.get("WEATHER_PROFILE_DIR")
    if weather_dir_str:
        weather_dir = _resolve_param_path(weather_dir_str, "weather directory")
        if not weather_dir.exists():
            missing_files.append(f"Weather directory: {weather_dir}")
    
    # Check rates file
    rates_file_str = params.get("RATES_FILE")
    if rates_file_str:
        rates_file = _resolve_param_path(rates_file_str, "rates file")
        if not rates_file.exists():
            missing_files.append(f"Rates file: {rates_file}")
    
    # Check planet files (if using planet executables)
    planet_dir_str = params.get("PLANET_DIR")
    planet_root_str = params.get("PLANET_ROOT")
    executable = params.get("EXECUTABLE", "")
    
    if planet_dir_str and planet_root_str and ("croin" in executable or "planet" in executable.lower()):
        planet_dir = _resolve_param_path(planet_dir_str, "planet directory")
        if not planet_dir.exists():
            missing_files.append(f"Planet directory: {planet_dir}")
    
    if missing_files:
        raise SmokeTestError(
            f"Missing input files:\n" +
            "\n".join(f"  - {f}" for f in missing_files) +
            "\n\nCheck that all file paths in your parameter file are correct."
        )


def verify_psf_files(params: Dict[str, str]) -> None:
    """Validate that detector PSF binaries are present and match detector expectations."""
    obs_dir_str = params.get("OBSERVATORY_DIR")
    obs_list_str = params.get("OBSERVATORY_LIST")
    if not obs_dir_str or not obs_list_str:
        return

    obs_dir = _resolve_param_path(obs_dir_str, "observatory directory")
    obs_list = _resolve_param_path(f"{obs_dir_str}/{obs_list_str}", "observatory list")
    if not obs_list.exists():
        raise SmokeTestError(f"Observatory list not found: {obs_list}")

    for raw_entry in obs_list.read_text(encoding="utf-8").splitlines():
        entry = raw_entry.strip()
        if not entry or entry.startswith("#"):
            continue

        observatory_path = (obs_dir / entry).resolve()
        if not observatory_path.is_file():
            raise SmokeTestError(
                f"Observatory configuration '{entry}' referenced in {obs_list.name} "
                f"was not found at {observatory_path}"
            )

        observatory_cfg = _parse_key_value_file(observatory_path)
        detector_rel = observatory_cfg.get("DETECTOR")
        if not detector_rel:
            # Some observatories embed detector parameters directly.
            continue

        detector_path = (obs_dir / detector_rel).resolve()
        if not detector_path.is_file():
            raise SmokeTestError(
                f"Detector file '{detector_rel}' referenced in {observatory_path.name} "
                f"was not found at {detector_path}"
            )

        detector_cfg = _parse_key_value_file(detector_path)
        psf_value = detector_cfg.get("PSFFILE", "").strip()
        if not psf_value:
            # Detector uses an analytic PSF; nothing to validate.
            continue

        psf_path = Path(psf_value)
        if not psf_path.is_absolute():
            psf_path = (REPO_ROOT / psf_value).resolve()
        if not psf_path.is_file():
            raise SmokeTestError(
                f"PSF file '{psf_value}' referenced in detector {detector_path.name} "
                f"was not found at {psf_path}"
            )

        try:
            psf_bytes = psf_path.read_bytes()
        except OSError as exc:
            raise SmokeTestError(f"Unable to read PSF file {psf_path}: {exc}") from exc

        try:
            version, ivars, dvars, header_bytes, total_size = _read_psf_header(psf_bytes)
        except struct.error as exc:
            raise SmokeTestError(f"PSF file {psf_path} is not a valid binary PSF: {exc}") from exc

        Nkern, nsub, psf_size, spsf_size, ipsf_size = ivars
        pixscale, _, _, psfh = dvars

        if version != 101:
            raise SmokeTestError(
                f"PSF file {psf_path} was produced with PSF class version {version}; "
                "expected version 101. Regenerate it with the current software."
            )

        expected_nkern = _safe_int(detector_cfg.get("KERNSIZE"))
        if expected_nkern is not None and Nkern != expected_nkern:
            raise SmokeTestError(
                f"PSF file {psf_path} reports Nkern={Nkern}, but detector {detector_path.name} "
                f"expects KERNSIZE={expected_nkern}. Regenerate the PSF with matching kernel size."
            )

        expected_nsub = _safe_int(detector_cfg.get("SUBPIX"))
        if expected_nsub is not None and nsub != expected_nsub:
            raise SmokeTestError(
                f"PSF file {psf_path} reports Nsub={nsub}, but detector {detector_path.name} "
                f"expects SUBPIX={expected_nsub}. Regenerate the PSF with matching subpixel sampling."
            )

        expected_pixscale = _safe_float(detector_cfg.get("PIXELSCALE"))
        if expected_pixscale is not None and not math.isclose(
            pixscale, expected_pixscale, rel_tol=1e-8, abs_tol=1e-10
        ):
            raise SmokeTestError(
                f"PSF file {psf_path} encodes PIXSCALE={pixscale}, but detector {detector_path.name} "
                f"expects PIXELSCALE={expected_pixscale}. Regenerate the PSF to match the detector."
            )

        if expected_pixscale is not None and expected_nsub is not None:
            expected_psfh = expected_pixscale / expected_nsub
            if not math.isclose(psfh, expected_psfh, rel_tol=1e-8, abs_tol=1e-10):
                raise SmokeTestError(
                    f"PSF file {psf_path} encodes step size {psfh}, but PIXELSCALE/SUBPIX = {expected_psfh}. "
                    "This usually indicates the PSF was generated with inconsistent inputs."
                )

        if expected_nkern is not None and expected_nsub is not None:
            expected_psf_size = (2 * expected_nkern + 1) ** 2 * (expected_nsub ** 2)
            if psf_size != expected_psf_size:
                raise SmokeTestError(
                    f"PSF file {psf_path} contains {psf_size} kernel samples, but detector {detector_path.name} "
                    f"requires {(2 * expected_nkern + 1)}×{(2 * expected_nkern + 1)} "
                    f"pixels with {expected_nsub}×{expected_nsub} subpixel sampling "
                    f"({expected_psf_size} samples). Regenerate the PSF."
                )

            small_nkern = expected_nkern // 2
            expected_spsf_size = (2 * small_nkern + 1) ** 2 * (expected_nsub ** 2)
            if spsf_size != expected_spsf_size:
                raise SmokeTestError(
                    f"PSF file {psf_path} has {spsf_size} samples in its small-kernel array, "
                    f"but detector {detector_path.name} expects {expected_spsf_size}."
                )

        expected_bytes = header_bytes + 8 * (psf_size + spsf_size + ipsf_size)
        if total_size != expected_bytes:
            raise SmokeTestError(
                f"PSF file {psf_path} has size {total_size:,} bytes, but its header implies "
                f"{expected_bytes:,} bytes. The file may be truncated or corrupted."
            )

        try:
            rel_path = psf_path.relative_to(REPO_ROOT).as_posix()
        except ValueError:
            rel_path = psf_path.as_posix()
        expected_hash = _CANONICAL_PSF_HASHES.get(rel_path)
        if expected_hash:
            digest = hashlib.sha256(psf_bytes).hexdigest()
            if digest != expected_hash:
                warnings.warn(
                    f"PSF file {psf_path} does not match the canonical SHA256 ({expected_hash}); "
                    "continuing because header metadata matches detector configuration."
                )


def verify_sequence_has_observations(params: Dict[str, str]) -> None:
    """Verify observing sequence files have at least one observation (Nstack > 0)."""
    obs_dir_str = params.get("OBSERVATORY_DIR")
    obs_list_str = params.get("OBSERVATORY_LIST")
    
    if not obs_dir_str or not obs_list_str:
        # Can't validate without observatory info
        return
    
    obs_dir = _resolve_param_path(obs_dir_str, "observatory directory")
    obs_list = _resolve_param_path(f"{obs_dir_str}/{obs_list_str}", "observatory list")
    
    if not obs_list.exists():
        # Observatory list doesn't exist, will fail elsewhere
        return
    
    # Read observatory list to find .observatory files
    lines = obs_list.read_text(encoding="utf-8").splitlines()
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        
        # Observatory list format is just filenames
        obs_file = obs_dir / line
        if not obs_file.exists():
            continue
        
        # Read observatory file to find sequence file
        obs_content = obs_file.read_text(encoding="utf-8")
        
        # Look for OBSERVATION_SEQUENCE keyword
        for obs_line in obs_content.splitlines():
            obs_line = obs_line.strip()
            if obs_line.startswith("OBSERVATION_SEQUENCE"):
                parts = obs_line.split()
                if len(parts) >= 2:
                    seq_file = obs_dir / parts[1]
                    
                    if not seq_file.exists():
                        raise SmokeTestError(
                            f"Sequence file not found: {seq_file}\n"
                            f"  Referenced in observatory file: {obs_file.name}"
                        )
                    
                    # Check sequence file has at least one observation
                    seq_lines = seq_file.read_text(encoding="utf-8").splitlines()
                    has_observation = False
                    
                    for seq_line in seq_lines:
                        seq_line = seq_line.strip()
                        if not seq_line or seq_line.startswith("#"):
                            continue
                        if seq_line.startswith("BEGIN_REPEAT") or seq_line.startswith("END_REPEAT"):
                            continue
                        
                        # Parse sequence line: Field Nstack Texp Sum Description
                        parts = seq_line.split()
                        if len(parts) >= 2:
                            try:
                                nstack = int(parts[1])
                                if nstack > 0:
                                    has_observation = True
                                    break
                            except ValueError:
                                continue
                    
                    if not has_observation:
                        raise SmokeTestError(
                            f"Sequence file {seq_file.name} has no observations (no lines with Nstack > 0).\n"
                            f"  At least one line must have positive Nstack to indicate this observatory takes images.\n"
                            f"  Format: Field Nstack(+ve) Texp Sum Description"
                        )
                    
                    # Only check first sequence file found
                    return


__all__ = [
    "verify_binary_source_columns",
    "verify_catalog_alignment", 
    "verify_catalog_columns",
    "verify_input_files_exist",
    "verify_psf_files",
    "verify_nfilters_matches_catalogs",
    "verify_outputs",
    "verify_rates_file",
    "verify_sequence_has_observations",
    "verify_source_lens_compatibility",
    "verify_weather_coverage",
]
