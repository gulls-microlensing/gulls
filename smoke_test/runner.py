"""Command-line entry point for the gulls smoke test."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .constants import (
    BUILD_BIN_DEFAULT,
    CASE_CHOICES,
    CASE_EXEC_MAP,
    CASE_LABELS,
    CASE_LOOKUP,
    CASES,
    REPO_ROOT,
)
from .errors import SmokeTestError
from .execution import run_command
from .metrics import gather_case_metrics
from .plotting import plot_lightcurves
from .prep import PreparedCase, prepare_cases
from .validation import (
    verify_binary_source_columns,
    verify_catalog_alignment,
    verify_catalog_columns,
    verify_input_files_exist,
    verify_psf_files,
    verify_nfilters_matches_catalogs,
    verify_outputs,
    verify_rates_file,
    verify_sequence_has_observations,
    verify_source_lens_compatibility,
    verify_weather_coverage,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-bin",
        type=Path,
        default=BUILD_BIN_DEFAULT,
        help="Directory containing the gulls executables (default: %(default)s)",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Skip cleaning existing smoke_test/output directories before running.",
    )
    parser.add_argument(
        "--cases",
        nargs="*",
        choices=list(CASE_CHOICES),
        help="Subset of cases to run (accepts case labels or executable names; default: all).",
    )
    parser.add_argument(
        "--instance",
        default="0",
        help="Instance identifier passed via -s (default: %(default)s).",
    )
    parser.add_argument(
        "--field",
        type=int,
        default=0,
        help="Field index passed via -f (default: %(default)s). Use --field -1 to let gulls auto-select.",
    )
    parser.add_argument(
        "--exec-timeout",
        type=float,
        default=180.0,
        help="Seconds to wait for each executable before aborting (<=0 disables).",
    )
    parser.add_argument(
        "--ci",
        action="store_true",
        help="Run CI-optimized subset of tests (faster, essential cases only).",
    )
    return parser.parse_args(argv)


def _resolve_case_selection(raw_choices: Sequence[str] | None, ci_mode: bool = False) -> Tuple[Tuple[str, str, str], ...]:
    if not raw_choices:
        if ci_mode:
            # CI subset: essential tests only (std, binary source validation)
            return (
                ("smoke_std", "gulls_std.x", "smoke_std.prm"),
                ("smoke_std_binary", "gulls_std.x", "smoke_std_binary.prm"),
                ("smoke_fish", "gullsFish.x", "smoke_fish.prm"),
                ("smoke_fish_binary", "gullsFish.x", "smoke_std_binary.prm"),
                ("smoke_croin", "gulls_croin.x", "smoke_croin.prm"),
                ("smoke_croin_binary", "gulls_croin.x", "smoke_croin_binary.prm"),
            )
        return CASES

    ordered: List[Tuple[str, str, str]] = []
    for choice in raw_choices:
        if choice in CASE_LABELS:
            ordered.append(CASE_LOOKUP[choice])
        else:
            ordered.extend(CASE_EXEC_MAP.get(choice, []))

    seen: set[str] = set()
    deduped: List[Tuple[str, str, str]] = []
    for case in ordered:
        if case[0] in seen:
            continue
        deduped.append(case)
        seen.add(case[0])
    return tuple(deduped)


def _prepare_environment(keep_output: bool, prepared_cases: Sequence[PreparedCase]) -> None:
    output_roots = {case.output_root for case in prepared_cases}
    if not keep_output:
        for root in output_roots:
            if root.exists():
                shutil.rmtree(root)
            root.mkdir(parents=True, exist_ok=True)
    else:
        for root in output_roots:
            root.mkdir(parents=True, exist_ok=True)


def _generate_psf_files(build_bin: Path) -> None:
    """Generate PSF files needed for smoke tests."""
    psf_dir = REPO_ROOT / "smoke_test" / "assets" / "observatories"
    psf_binary_file = psf_dir / "WFI_PSF.psf"
    
    # Check if we already have a valid PSF file (should be ~68MB for subpixel sampling)
    if psf_binary_file.exists() and psf_binary_file.stat().st_size > 10_000_000:  # > 10MB
        print(f"Using existing PSF file: {psf_binary_file} ({psf_binary_file.stat().st_size:,} bytes)")
        return
    
    # Generate PSF using PSF class with Moffat function
    generate_moffat_psf = build_bin / "generateMoffatPSF"
    
    if not generate_moffat_psf.exists():
        raise SmokeTestError(f"PSF generator not found: {generate_moffat_psf}")
    
    # Generate binary PSF with subpixel sampling
    print("Generating PSF with subpixel sampling using PSF class...")
    cmd = [
        str(generate_moffat_psf),
        "0.2",    # fwhm (arcsec) - matches PSFFWHM
        "0.11",   # pixel_scale (arcsec) - matches PIXELSCALE
        str(psf_binary_file)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise SmokeTestError(f"PSF generation failed: {result.stderr}")
    
    if not psf_binary_file.exists():
        raise SmokeTestError(f"PSF file was not created: {psf_binary_file}")
    
    print(f"Generated PSF file: {psf_binary_file} ({psf_binary_file.stat().st_size:,} bytes)")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    build_bin = args.build_bin.resolve()
    if not build_bin.is_dir():
        raise SmokeTestError(f"Build directory not found: {build_bin}")

    selected = _resolve_case_selection(args.cases, args.ci)
    if not selected:
        print("No cases selected", file=sys.stderr)
        return 1

    espl_table = REPO_ROOT / "src" / "ESPL.tbl"
    if not espl_table.is_file():
        raise SmokeTestError(f"Missing ESPL.tbl at {espl_table}; copy it before running the smoke test.")

    # Generate PSF files needed for smoke tests
    _generate_psf_files(build_bin)

    prepared_cases, prep_failures = prepare_cases(build_bin, selected)
    if prep_failures:
        print("Smoke test setup issues:")
        for message in prep_failures:
            print(f" - {message}")

    if not prepared_cases:
        return 1
    
    # Validate catalogs: columns, compatibility, and binary-specific requirements
    for case in prepared_cases:
        try:
            verify_input_files_exist(case.params)
            verify_psf_files(case.params)
            verify_catalog_columns(case.params)
            verify_source_lens_compatibility(case.params)
            verify_binary_source_columns(case.params)
            verify_nfilters_matches_catalogs(case.params)
            verify_weather_coverage(case.params)
            verify_rates_file(case.params)
            verify_sequence_has_observations(case.params)
        except SmokeTestError as err:
            print(f"Validation failed for {case.label}:")
            print(f" - {err}")
            return 1

    env = os.environ.copy()
    base_dir = REPO_ROOT.as_posix() + "/"
    env["GULLS_BASE_DIR"] = base_dir
    env.setdefault("GULLS_STARS_DIR", base_dir)

    _prepare_environment(args.keep_output, prepared_cases)

    failures: List[str] = []

    for case in prepared_cases:
        if case.output_dir.exists() and not args.keep_output:
            shutil.rmtree(case.output_dir)
        case.output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [str(case.exe_path), "-i", str(case.exec_param_path), "-s", args.instance]
        if args.field is not None:
            cmd.extend(["-f", str(args.field)])

        print(f"\n=== Running {case.exe_name} ({case.label}) with {case.param_path.name} ===")
        result = run_command(cmd, env, args.exec_timeout, cwd=REPO_ROOT)

        print(result.stdout)
        if result.returncode != 0:
            failures.append(f"{case.label} exited with {result.returncode}")
            continue

        out_files = verify_outputs(case.output_dir)
        verify_catalog_alignment(out_files, case.params)
        summaries = gather_case_metrics(out_files)
        plot_lightcurves(case.output_dir, summaries, case.params)

    if failures:
        print("\nSmoke test failed:")
        for message in failures:
            print(f" - {message}")
        return 1

    print("\nAll smoke test cases completed successfully.")
    return 0


__all__ = ["main", "parse_args"]
