#!/usr/bin/env python3
"""
Validate observatory configuration files for required parameters.
"""

import sys
import os

def validate_observatory_file(filepath):
    """Validate an observatory file has all required parameters."""
    
    required_params = [
        "NAME", "LATITUDE", "LONGITUDE", "ALTITUDE", "READ_OHEAD", "NFIELDS",
        "WEATHER_PROFILE", "FIELDCENTRES", "NPIX_X", "NPIX_Y", "PIXELSIZE",
        "PRIMARY", "BLOCKAGE", "SPACE", "FILTER", "OBSERVATION_SEQUENCE",
        "DETECTOR", "THROUGHPUT", "REFERENCE_TEXP", "REFERENCE_NSTACK",
        "ORBIT", "PHOTOMETRY", "EXTCOEFF", "SKY_BACKGROUND"
    ]
    
    optional_params = [
        "MOON_AVOID", "ALT_LIMIT"
    ]
    
    if not os.path.exists(filepath):
        print(f"ERROR: Observatory file not found: {filepath}")
        return False
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    missing_required = []
    missing_optional = []
    
    for param in required_params:
        if param not in content:
            missing_required.append(param)
    
    for param in optional_params:
        if param not in content:
            missing_optional.append(param)
    
    if missing_required:
        print(f"ERROR: Missing required parameters in {filepath}:")
        for param in missing_required:
            print(f"  - {param}")
        return False
    
    if missing_optional:
        print(f"WARNING: Missing optional parameters in {filepath}:")
        for param in missing_optional:
            print(f"  - {param}")
        print("These will use default values.")
    
    print(f"✓ Observatory file {filepath} is valid")
    return True

def main():
    if len(sys.argv) != 2:
        print("Usage: validate_observatory.py <observatory_file>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    if validate_observatory_file(filepath):
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
