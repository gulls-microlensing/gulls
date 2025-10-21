#!/usr/bin/env python3
"""Extract all Gulls parameters from readParamfile.cpp for documentation."""

import re
import sys
from pathlib import Path

def extract_parameters():
    """Extract parameter definitions from readParamfile.cpp."""
    paramfile_path = Path("src/readParamfile.cpp")
    
    if not paramfile_path.exists():
        print(f"Error: {paramfile_path} not found")
        return
    
    with open(paramfile_path, 'r') as f:
        content = f.read()
    
    # Extract default values from the defaults map
    defaults_pattern = r'{"([^"]+)","([^"]*)"}'
    defaults = dict(re.findall(defaults_pattern, content))
    
    # Extract parameter assignments
    param_pattern = r'Paramfile->(\w+)\s*=\s*[^;]+;'
    params = re.findall(param_pattern, content)
    
    # Extract stod/stoi calls for type information
    type_pattern = r'Paramfile->(\w+)\s*=\s*(stod|stoi|stof)\([^)]+\)'
    types = dict(re.findall(type_pattern, content))
    
    print("# Gulls Parameter Reference")
    print()
    print("This document lists all parameters available in Gulls parameter files.")
    print()
    
    # Group parameters by category
    categories = {
        "Basic Configuration": [
            "RUN_NAME", "OUTPUT_DIR", "FINAL_DIR", "EXECUTABLE"
        ],
        "Input Directories": [
            "OBSERVATORY_DIR", "OBSERVATORY_LIST", "WEATHER_PROFILE_DIR",
            "STARFIELD_DIR", "STARFIELD_LIST", "SOURCE_DIR", "SOURCE_LIST",
            "LENS_DIR", "LENS_LIST", "PLANET_DIR", "PLANET_ROOT", "RATES_FILE"
        ],
        "Output Control": [
            "OUTPUT_LC", "OUTPUT_IMAGES", "PRETTY_PICS", "PRETTY_PICS_DIMENSIONS",
            "OUTPUT_ONERR", "OUTPUT_ONDET", "OUTPUT_ONALL", "VERBOSITY"
        ],
        "Simulation Parameters": [
            "SET_RANDOM_SEED_TO_CLOCK", "RANDOM_SEED", "SIMULATION_ZERO_TIME",
            "SUBRUNSIZE", "NUM_SIM_DAYS", "REPEAT_SEQUENCE"
        ],
        "Physical Parameters": [
            "NFILTERS", "AMIN", "LARGEPSFMAG", "MIN_CHISQUARED", "U0MAX",
            "PARALLAX", "LENS_LIGHT", "ASTROMETRY_ON", "ASTROMETRIC_SYS_FLOOR"
        ],
        "Lightcurve Generation": [
            "LC_GEN", "LD_GAMMA", "VBM_RELTOL", "VBM_ABSTOL", "LC_TIMEOUT"
        ],
        "Multiplicity": [
            "MULTIPLE_SOURCES", "MULTIPLE_LENSES"
        ],
        "Observatory Groups": [
            "OBS_GROUPS", "OBS_GROUP_NAMES", "ERROR_SCALING"
        ]
    }
    
    for category, param_list in categories.items():
        print(f"## {category}")
        print()
        
        for param in param_list:
            if param in defaults:
                default_val = defaults[param]
                param_type = types.get(param.lower(), "string")
                
                print(f"### {param}")
                print(f"- **Type**: {param_type}")
                print(f"- **Default**: `{default_val}`")
                print(f"- **Description**: [TODO - needs manual documentation]")
                print()
    
    # List any parameters not in categories
    all_categorized = set()
    for param_list in categories.values():
        all_categorized.update(param_list)
    
    uncategorized = set(defaults.keys()) - all_categorized
    if uncategorized:
        print("## Other Parameters")
        print()
        for param in sorted(uncategorized):
            default_val = defaults[param]
            param_type = types.get(param.lower(), "string")
            print(f"### {param}")
            print(f"- **Type**: {param_type}")
            print(f"- **Default**: `{default_val}`")
            print(f"- **Description**: [TODO - needs manual documentation]")
            print()

if __name__ == "__main__":
    extract_parameters()
