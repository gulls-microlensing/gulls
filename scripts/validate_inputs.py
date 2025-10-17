#!/usr/bin/env python3
"""Validate GULLS input files before running simulations."""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Any

def validate_sources(filepath: Path) -> List[str]:
    """Validate source catalog format."""
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
        required_cols = ['RA2000.0', 'DEC2000.0', 'Dist', 'Mass', 'Radius']
        
        for col in required_cols:
            if col not in header:
                errors.append(f"Missing required column: {col}")
        
        # Check for data
        if len(lines) < 2:
            errors.append("No data rows found")
        
        # Basic data validation
        for i, line in enumerate(lines[1:], 2):
            if line.strip().startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) != len(header):
                errors.append(f"Row {i}: Expected {len(header)} columns, got {len(parts)}")
                continue
            
            # Check numeric values
            try:
                ra = float(parts[header.index('RA2000.0')])
                dec = float(parts[header.index('DEC2000.0')])
                dist = float(parts[header.index('Dist')])
                mass = float(parts[header.index('Mass')])
                radius = float(parts[header.index('Radius')])
                
                if not (0 <= ra <= 360):
                    errors.append(f"Row {i}: RA must be 0-360 degrees")
                if not (-90 <= dec <= 90):
                    errors.append(f"Row {i}: DEC must be -90 to 90 degrees")
                if dist <= 0:
                    errors.append(f"Row {i}: Distance must be positive")
                if mass <= 0:
                    errors.append(f"Row {i}: Mass must be positive")
                if radius <= 0:
                    errors.append(f"Row {i}: Radius must be positive")
                    
            except (ValueError, IndexError) as e:
                errors.append(f"Row {i}: Invalid numeric data: {e}")
    
    except Exception as e:
        errors.append(f"Error reading source file: {e}")
    
    return errors

def validate_lenses(filepath: Path) -> List[str]:
    """Validate lens catalog format."""
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
        required_cols = ['RA2000.0', 'DEC2000.0', 'Dist', 'Mass', 'mul', 'mub']
        
        for col in required_cols:
            if col not in header:
                errors.append(f"Missing required column: {col}")
        
        # Check for data
        if len(lines) < 2:
            errors.append("No data rows found")
        
        # Basic data validation
        for i, line in enumerate(lines[1:], 2):
            if line.strip().startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) != len(header):
                errors.append(f"Row {i}: Expected {len(header)} columns, got {len(parts)}")
                continue
            
            # Check numeric values
            try:
                ra = float(parts[header.index('RA2000.0')])
                dec = float(parts[header.index('DEC2000.0')])
                dist = float(parts[header.index('Dist')])
                mass = float(parts[header.index('Mass')])
                
                if not (0 <= ra <= 360):
                    errors.append(f"Row {i}: RA must be 0-360 degrees")
                if not (-90 <= dec <= 90):
                    errors.append(f"Row {i}: DEC must be -90 to 90 degrees")
                if dist <= 0:
                    errors.append(f"Row {i}: Distance must be positive")
                if mass <= 0:
                    errors.append(f"Row {i}: Mass must be positive")
                    
            except (ValueError, IndexError) as e:
                errors.append(f"Row {i}: Invalid numeric data: {e}")
    
    except Exception as e:
        errors.append(f"Error reading lens file: {e}")
    
    return errors

def validate_parameter_file(filepath: Path) -> List[str]:
    """Validate parameter file format."""
    errors = []
    
    if not filepath.exists():
        return [f"Parameter file not found: {filepath}"]
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        required_params = [
            'RUN_NAME', 'OUTPUT_DIR', 'EXECUTABLE',
            'OBSERVATORY_DIR', 'OBSERVATORY_LIST',
            'SOURCE_DIR', 'SOURCE_LIST',
            'LENS_DIR', 'LENS_LIST'
        ]
        
        found_params = set()
        
        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if '=' not in line:
                errors.append(f"Line {i}: Invalid format (missing '=')")
                continue
            
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip()
            
            found_params.add(key)
            
            if not value:
                errors.append(f"Line {i}: Parameter {key} has no value")
        
        for param in required_params:
            if param not in found_params:
                errors.append(f"Missing required parameter: {param}")
    
    except Exception as e:
        errors.append(f"Error reading parameter file: {e}")
    
    return errors

def main():
    parser = argparse.ArgumentParser(description="Validate GULLS input files")
    parser.add_argument("--sources", type=Path, help="Source catalog file")
    parser.add_argument("--lenses", type=Path, help="Lens catalog file")
    parser.add_argument("--params", type=Path, help="Parameter file")
    parser.add_argument("--all", action="store_true", help="Validate all files in current directory")
    
    args = parser.parse_args()
    
    all_errors = []
    
    if args.all:
        # Find files automatically
        sources_files = list(Path(".").glob("**/*.sources"))
        lens_files = list(Path(".").glob("**/*.lenses"))
        param_files = list(Path(".").glob("**/*.prm"))
        
        for sources_file in sources_files:
            all_errors.extend(validate_sources(sources_file))
        
        for lens_file in lens_files:
            all_errors.extend(validate_lenses(lens_file))
        
        for param_file in param_files:
            all_errors.extend(validate_parameter_file(param_file))
    
    else:
        if args.sources:
            all_errors.extend(validate_sources(args.sources))
        
        if args.lenses:
            all_errors.extend(validate_lenses(args.lenses))
        
        if args.params:
            all_errors.extend(validate_parameter_file(args.params))
    
    if all_errors:
        print("Validation errors found:")
        for error in all_errors:
            print(f"  - {error}")
        return 1
    else:
        print("All files validated successfully!")
        return 0

if __name__ == "__main__":
    sys.exit(main())
