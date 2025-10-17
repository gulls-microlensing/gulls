# Gulls Smoke Test Suite

This directory contains a minimal end-to-end smoke test for the three gulls executables: `gulls_std`, `gulls_croin`, and `gullsFish`. The smoke test verifies that the executables can run successfully with synthetic input data and produce expected outputs including lightcurves, reports, and visualizations.

## Overview

The smoke test uses synthetic catalogs and simplified observatory configurations to quickly validate that:
- All three executables compile and run without errors
- Lightcurve generation completes within timeout limits
- Photometry and astrometry outputs are generated correctly
- Visualization plots can be created from the output data

## Directory Structure

```
smoke_test/
├── README.md                  # This file
├── run_smoke_test.py          # Main test runner script
├── parameterfiles/            # Parameter files for each executable
│   ├── smoke_std.prm
│   ├── smoke_croin.prm
│   └── smoke_fish.prm
├── assets/                    # Test input data
│   ├── lenses/                # Synthetic lens catalogs
│   ├── sources/               # Synthetic source catalogs
│   ├── starfields/            # Starfield definitions
│   ├── observatories/         # Observatory configurations
│   ├── rates/                 # Event rate files
│   ├── weather/               # Weather profiles
│   └── planets/               # Planet configurations
└── output/                    # Test outputs (generated at runtime)
    ├── std/
    ├── croin/
    └── fish/
```

## Catalog Helper Tool

Use `dat_tool.py` to inspect and adjust the whitespace-delimited catalogs under `assets/`.

- Inspect column layout with optional sample values:
  ```bash
  python dat_tool.py show assets/lenses/smoke_lens_catalog.dat
  ```
- Scale or offset a numeric column (writes to a new file unless `--inplace` is provided):
  ```bash
  python dat_tool.py apply assets/lenses/smoke_lens_catalog.dat \
      --column Mass --operation multiply --value 1.5 \
      --output assets/lenses/smoke_lens_catalog_scaled.dat
  ```
- Use `--format` to control numeric formatting (default `"{value:.7e}"`) and `--preview`
  to choose how many modified rows are echoed.

## Requirements

- Built gulls executables in `build/bin/` (see main README for build instructions)
- Python 3.7+ with the following packages:
  - `numpy`
  - `pandas`
  - `matplotlib` (optional, for plotting)

### Setting Up Python Dependencies

#### Using Conda (Recommended)
A conda environment file is provided for easy setup:
```bash
conda env create -f smoke_test/smoke.yml
conda activate smoke
```

#### Using pip
Alternatively, install packages with pip:
```bash
pip install numpy pandas matplotlib
```

## Running the Smoke Test

### Basic Usage

Run all three test cases:
```bash
python smoke_test/run_smoke_test.py
```

Run a specific executable:
```bash
python smoke_test/run_smoke_test.py --cases gulls_std
```

Run with custom timeout (default: 180 seconds):
```bash
python smoke_test/run_smoke_test.py --exec-timeout 120
```

### Command-Line Options

- `--build-bin PATH`: Directory containing executables (default: `build/bin`)
- `--keep-output`: Skip cleaning existing output directories before running
- `--cases CASE [CASE ...]`: Subset of runs to execute. Accepts either executable names (`gulls_std`, `gulls_croin`, `gullsFish`) or case labels (`std-single`, `std-binary`, `std-heavy`, `croin-single`, `croin-binary`, `croin-heavy`, `fish-single`, `fish-binary`, `fish-heavy`). Case-specific run names are appended automatically (for example, `smoke_std_std-heavy`) so the heavy scenarios do not overwrite the baseline outputs.
- `--instance ID`: Instance identifier passed via `-s` flag (default: `0`)
- `--field N`: Field index passed via `-f` flag (default: `0`; use `-1` for auto-select)
- `--exec-timeout SECONDS`: Timeout per executable (default: `180`; `<=0` disables)

### Examples

```bash
# Test only gulls_std with 2-minute timeout
python smoke_test/run_smoke_test.py --cases gulls_std --exec-timeout 120

# Test all executables and keep previous outputs
python smoke_test/run_smoke_test.py --keep-output

# Test with custom build directory
python smoke_test/run_smoke_test.py --build-bin /path/to/custom/build/bin
```

## Test Configuration

### Simulation Parameters

The smoke test uses simplified parameters for quick execution:
- **Simulation length**: 100 days
- **Lightcurve timeout**: 10 seconds (enforced on VBMicrolensing)
- **Single observatory**: "SmokeScope" with F184 filter
- **Exposure time**: 5 seconds
- **Astrometry**: Enabled with 0.1 mas systematic floor
- **6 synthetic sources and 6 synthetic lenses** per test

### Observatory Configuration

The test observatory (`smoke.observatory`) is configured with:
- **Filter**: F184 (index 6, ~2 μm)
- **Detector**: 10M electron full well capacity
- **Exposure**: 5 seconds per observation
- **Zero magnitude**: 20.0
- **Sequence**: 100 days of observations

## Expected Outputs

For each test case, the following outputs are generated:

### 1. Summary Reports (`.out`)
Located in `smoke_test/output/{std,croin,fish}/smoke_{std,croin,fish}/`

Contains run statistics, parameter values, and execution timing.

### 2. Lightcurve Files (`.lc`)
ASCII tables with columns including:
- `Simulation_time`: Observation epoch (days)
- `measured_relative_flux`: Photometric flux (relative to baseline)
- `measured_relative_flux_error`: Photometric error
- `true_relative_flux`: True simulated flux
- Astrometric centroids in pixel coordinates
- Astrometric shifts in milliarcseconds (N/E)
- Sky positions in RA/Dec (degrees)
- Astrometric errors
- Parallax information
- Source and lens positions

### 3. Visualization Plots (`.png`)
If matplotlib is available, the script generates diagnostic plots showing:
- **Top-left**: Lightcurve (flux vs time) with measured points, true flux line, and baseline
- **Top-right**: Astrometric position in RA/Dec with error bars
- **Bottom-left**: Astrometric trajectory in North/East coordinates
- **Bottom-right**: Astrometric centroid shifts vs time

## Validation Criteria

The smoke test passes when:
1. All selected executables complete without errors (exit code 0)
2. At least one `.out` summary file is created and non-empty
3. At least one `.lc` lightcurve file is generated
4. No timeouts occur during execution

## Troubleshooting

### "Missing executable" error
- Ensure you've built the project: `cmake --build build`
- Check that executables exist in `build/bin/`

### "Missing ESPL.tbl" error
- Copy the VBMicrolensing ESPL table: `cp VBMicrolensing/ESPL.tbl src/`

### Timeout errors
- Increase timeout: `--exec-timeout 300`
- Check for infinite loops in lightcurve generation
- Verify VBMicrolensing `LC_TIMEOUT` is enforced (see `src/pllxLightcurveGenerator.cpp`)

### Plotting errors
- Install required packages: `pip install numpy pandas matplotlib`
- Run without plotting by uninstalling matplotlib temporarily

### "Negative error values" warning
- This was fixed by using pandas with named columns
- Verify column names match between code and output files

## Technical Details

### VBMicrolensing Timeout Enforcement
The smoke test relies on timeout enforcement in VBMicrolensing to prevent infinite hangs during lightcurve generation. The timeout is set via `LC_TIMEOUT=10.0` in parameter files and enforced in `src/pllxLightcurveGenerator.cpp`.

### Pandas-Based Column Access
The plotting function uses pandas DataFrames with named columns for robust data access. This prevents column index misalignment issues and makes the code more maintainable.

### Data Format
Lightcurve files use whitespace-delimited ASCII format with:
- Comment lines starting with `#`
- Header line with column names
- Data rows with numerical values

## Performance

Typical execution times (MacBook Pro, M1):
- **gulls_std**: ~10 seconds
- **gulls_croin**: ~10 seconds  
- **gullsFish**: ~10 seconds
- **Total** (all three): ~30 seconds

## Contributing

When modifying the smoke test:
1. Keep simulation parameters minimal for fast execution
2. Ensure all three executables are tested
3. Verify plots are generated correctly
4. Update this README if you add new features or change behavior
