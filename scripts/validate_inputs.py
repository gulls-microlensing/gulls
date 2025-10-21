#!/usr/bin/env python3
"""
Validate Gulls input files before running simulations.

This script checks:
- Parameter file format and required parameters
- Source and lens catalog column headers
- Valid source/lens distance pairs
- Reasonable numerical ranges
- Binary source requirements (when MULTIPLE_SOURCES=1)

Usage:
    python validate_inputs.py <parameter_file.prm>
    python validate_inputs.py --sources <file> --lenses <file>
    python validate_inputs.py --all

Examples:
    # Validate using parameter file (recommended - most comprehensive)
    python validate_inputs.py my_simulation.prm
    
    # Validate specific files
    python validate_inputs.py --sources my_sources.dat --lenses my_lenses.dat
    
    # Auto-discover and validate all input files in current directory
    python validate_inputs.py --all
"""

import argparse
import sys
from pathlib import Path
from typing import List

# Add repo root to path to import smoke_test as a package
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

try:
    from smoke_test.validation import (
        verify_binary_source_columns,
        verify_catalog_columns,
        verify_input_files_exist,
        verify_nfilters_matches_catalogs,
        verify_rates_file,
        verify_sequence_has_observations,
        verify_source_lens_compatibility,
        verify_weather_coverage,
    )
    from smoke_test.errors import SmokeTestError
    USE_SMOKE_TEST_VALIDATION = True
except ImportError as e:
    USE_SMOKE_TEST_VALIDATION = False
    print(f"Warning: Could not import smoke test validation ({e}). Using basic checks only.")


def read_parameter_file(param_file: Path) -> dict:
    """Parse a gulls parameter file into a dictionary."""
    params = {}
    
    with param_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue
            
            # Split on = sign
            if "=" in line:
                key, value = line.split("=", 1)
                params[key.strip()] = value.strip()
    
    return params


def validate_parameter_file(filepath: Path) -> List[str]:
    """Validate parameter file format and required parameters."""
    errors = []
    
    if not filepath.exists():
        return [f"Parameter file not found: {filepath}"]
    
    try:
        params = read_parameter_file(filepath)
        
        required_params = [
            'RUN_NAME', 'OUTPUT_DIR', 'EXECUTABLE',
            'OBSERVATORY_DIR', 'OBSERVATORY_LIST',
            'SOURCE_DIR', 'SOURCE_LIST',
            'LENS_DIR', 'LENS_LIST'
        ]
        
        for param in required_params:
            if param not in params:
                errors.append(f"Missing required parameter: {param}")
            elif not params[param]:
                errors.append(f"Parameter {param} has no value")
    
    except Exception as e:
        errors.append(f"Error reading parameter file: {e}")
    
    return errors


def validate_sources_basic(filepath: Path) -> List[str]:
    """Basic validation of source catalog format (fallback when smoke test unavailable)."""
    errors = []
    
    if not filepath.exists():
        return [f"Source file not found: {filepath}"]
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            return ["Source file is empty"]
        
        # Check header
        header = lines[0].strip().split()
        required_cols = {'mul', 'mub', 'Mass', 'Radius', 'Dist'}
        missing = required_cols - set(header)
        
        if missing:
            errors.append(f"Missing required columns: {', '.join(sorted(missing))}")
        
        # Check for data
        if len(lines) < 2:
            errors.append("No data rows found")
    
    except Exception as e:
        errors.append(f"Error reading source file: {e}")
    
    return errors


def validate_lenses_basic(filepath: Path) -> List[str]:
    """Basic validation of lens catalog format (fallback when smoke test unavailable)."""
    errors = []
    
    if not filepath.exists():
        return [f"Lens file not found: {filepath}"]
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            return ["Lens file is empty"]
        
        # Check header
        header = lines[0].strip().split()
        required_cols = {'mul', 'mub', 'Mass', 'Dist'}
        missing = required_cols - set(header)
        
        if missing:
            errors.append(f"Missing required columns: {', '.join(sorted(missing))}")
        
        # Check for data
        if len(lines) < 2:
            errors.append("No data rows found")
    
    except Exception as e:
        errors.append(f"Error reading lens file: {e}")
    
    return errors


def check_numerical_libraries() -> List[str]:
    """Check if GSL fallbacks are being used instead of Numerical Recipes."""
    warnings = []
    
    # Check if random.cpp exists and contains GSL fallback code
    random_cpp = REPO_ROOT / "src" / "classes" / "random.cpp"
    if random_cpp.exists():
        content = random_cpp.read_text(encoding="utf-8")
        if "gsl_rng" in content and "Fallback implementation" in content:
            warnings.append(
                "⚠ Warning: Using GSL fallback implementations for random number generation.\n"
                "  The actual Numerical Recipes implementations are preferred for production runs.\n"
                "  GSL fallbacks are provided for CI/testing purposes only."
            )
    
    # Check if zroots2.cpp exists and contains GSL fallback code
    zroots_cpp = REPO_ROOT / "src" / "classes" / "zroots2.cpp"
    if zroots_cpp.exists():
        content = zroots_cpp.read_text(encoding="utf-8")
        if "gsl_poly_complex_solve" in content and "Fallback implementation" in content:
            warnings.append(
                "⚠ Warning: Using GSL fallback implementations for polynomial root finding.\n"
                "  The actual Numerical Recipes implementations are preferred for production runs.\n"
                "  GSL fallbacks are provided for CI/testing purposes only."
            )
    
    return warnings


def validate_with_smoke_test(param_file: Path) -> List[str]:
    """Run comprehensive validation using smoke test validation functions."""
    errors = []
    
    try:
        params = read_parameter_file(param_file)
        
        # Check for GSL fallbacks first
        warnings = check_numerical_libraries()
        if warnings:
            print("\n" + "=" * 70)
            for warning in warnings:
                print(warning)
            print("=" * 70)
        
            print("\n1. Checking parameter file format...")
            param_errors = validate_parameter_file(param_file)
            if param_errors:
                errors.extend(param_errors)
                return errors  # Can't proceed without valid params
            print("   ✓ Parameter file format valid")
            
            print("\n2. Checking input files exist...")
            verify_input_files_exist(params)
            print("   ✓ All input files found")
            
            print("\n3. Checking required catalog columns...")
            verify_catalog_columns(params)
            print("   ✓ All required columns present")
            
            print("\n4. Checking source/lens compatibility...")
            verify_source_lens_compatibility(params)
            print("   ✓ Valid source/lens pairs exist")
            print("   ✓ Distance values are reasonable")
            
            print("\n5. Checking binary source requirements...")
            multiple_sources = params.get("MULTIPLE_SOURCES", "0").strip()
            if multiple_sources in ("1", "1.0"):
                verify_binary_source_columns(params)
                print("   ✓ Binary source columns present")
            else:
                print("   ⊘ Skipped (MULTIPLE_SOURCES not enabled)")
            
            print("\n6. Checking NFILTERS matches catalog format...")
            verify_nfilters_matches_catalogs(params)
            print("   ✓ NFILTERS matches magnitude columns")
            
            print("\n7. Checking weather file coverage...")
            verify_weather_coverage(params)
            print("   ✓ Weather file covers simulation duration")
            
            print("\n8. Checking rates file validity...")
            verify_rates_file(params)
            print("   ✓ Rates file parameters are valid")
            
            print("\n9. Checking observing sequence has observations...")
            verify_sequence_has_observations(params)
            print("   ✓ Sequence file contains observations")
        
    except SmokeTestError as e:
        errors.append(str(e))
    except Exception as e:
        errors.append(f"Validation error: {e}")
    
    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate Gulls input files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("param_file", nargs="?", type=Path, 
                       help="Parameter file (.prm) - validates all catalogs referenced in it")
    parser.add_argument("--sources", type=Path, help="Source catalog file (basic validation only)")
    parser.add_argument("--lenses", type=Path, help="Lens catalog file (basic validation only)")
    parser.add_argument("--all", action="store_true", 
                       help="Auto-discover and validate all .prm files in current directory")
    
    args = parser.parse_args()
    
    all_errors = []
    
    if args.param_file and USE_SMOKE_TEST_VALIDATION:
        # Comprehensive validation using parameter file
        print(f"Validating catalogs specified in: {args.param_file}")
        print("=" * 70)
        all_errors = validate_with_smoke_test(args.param_file)
        
    elif args.all:
        # Find and validate all parameter files
        param_files = list(Path(".").glob("**/*.prm"))
        
        if not param_files:
            print("No .prm files found in current directory")
            return 1
        
        for param_file in param_files:
            print(f"\nValidating {param_file}...")
            print("=" * 70)
            
            if USE_SMOKE_TEST_VALIDATION:
                errors = validate_with_smoke_test(param_file)
            else:
                errors = validate_parameter_file(param_file)
            
            if errors:
                all_errors.extend([f"{param_file}: {e}" for e in errors])
    
    else:
        # Basic individual file validation (legacy mode)
        if args.sources:
            print(f"Validating source catalog: {args.sources}")
            all_errors.extend(validate_sources_basic(args.sources))
        
        if args.lenses:
            print(f"Validating lens catalog: {args.lenses}")
            all_errors.extend(validate_lenses_basic(args.lenses))
        
        if not args.sources and not args.lenses and not args.param_file:
            parser.print_help()
            return 1
    
    # Print results
    if all_errors:
        print("\n" + "=" * 70)
        print("✗ Validation failed:")
        for error in all_errors:
            print(f"  {error}")
        print("\n" + "=" * 70)
        print("Please fix the issues above and try again.")
        return 1
    else:
        print("\n" + "=" * 70)
        print("✓ All validations passed!")
        print("\nYour input files appear to be correctly formatted.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
