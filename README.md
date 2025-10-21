# Gulls: Gravitational Microlensing Simulation Suite

Gulls is a comprehensive simulation framework for modeling gravitational microlensing events with realistic observing conditions.

## Quick Start

### For Existing Users
**Everything you know still works!** Your existing workflow, parameter files, and scripts are unchanged.

### For New Users
Choose your preferred build method:

**Option 1: CMake (recommended)**
```bash
git clone <repository-url>
cd gulls
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

**Option 2: Traditional Makefile**
```bash
git clone <repository-url>
cd gulls
./configure.sh  # Required: sets up paths and environment
make clean
make
```

**Both methods produce identical executables!**

### Legacy Makefile Workflow

The historical Makefiles remain in `src/` for anyone who relies on the old process. They still expect manual setup (absolute paths, environment variables, etc.) and are **no longer recommended** for new users.

> **Note**: Some older Makefile targets reference Fortran sources that are not present in the repository. Only the `gulls_std`, `gulls_croin`, and `gullsFish` executables currently build successfully.

### Prerequisites
- **C++17** compatible compiler (GCC, Clang, or MSVC)
- **Fortran** compiler (gfortran recommended)
- **CMake** 3.10+ (for CMake build) or **Make** (for traditional build)
- **Python 3.7+** (for validation and testing - optional)

### Run a Simulation
```bash
# Your existing workflow (unchanged)
./bin/gulls_std.x -i parameter_file.prm -s 0 -f 0

# Or with CMake build
./build/bin/gulls_std.x -i parameter_file.prm -s 0 -f 0
```

### Optional: Validate Your Input Catalogs
```bash
# Check your input files before running (catches common errors)
python scripts/validate_inputs.py parameter_file.prm
```

### Optional: Test the Installation
```bash
# Run automated tests to verify everything works
python smoke_test/run_smoke_test.py --ci
```

## 🔧 Diagnostic Runs and Troubleshooting

### Diagnostic Configuration

For initial testing and debugging, use these settings in your parameter file:

```bash
PRETTY_PICS=1
PRETTY_PICS_DIMENSIONS=128  # Or larger for better diagnostics
OUTPUT_LC=1                 # Output every lightcurve
OUTPUT_IMAGES=1             # Save baseline and peak images
OUTPUT_ONALL=1              # Output all events (not just detections)
```

This produces maximum diagnostic output:
- **Lightcurves** (`.lc` files): `.det.lc` shows detectable events, `.all.lc` shows all events
- **Images**: One at peak magnification, one at baseline for each filter
- **Debug information**: Use `-d` flags for increasing verbosity

### Troubleshooting Common Issues

**White square images in ds9:**
1. Select "zscale" scaling in ds9
2. If still white, check for bad magnitudes in star/source/lens catalogs
3. Try larger `PRETTY_PICS_DIMENSIONS` (e.g., 512) to reduce impact of bright stars
4. Inspect color-magnitude diagrams of catalogs for anomalies

**"No valid fields were loaded":**
- Check that `NFILTERS` matches your observatory configuration
- Verify observatory files are in correct paths

**Simulation hangs or times out:**
- Check `LC_TIMEOUT` parameter (default 10 seconds may be too short for complex events)
- Increase for challenging binary configurations (e.g., `LC_TIMEOUT=600`)
- Review validation warnings about catalog issues

**GSL fallback warnings:**
- For production science, replace `src/classes/random.cpp` and `src/classes/zroots2.cpp` with licensed Numerical Recipes implementations
- GSL fallbacks are for CI/testing only
- **Binary releases include GSL fallbacks and are NOT suitable for production science**

### Image Analysis

Inspect images in ds9 using zscale to check for:
- **Realistic starfields** - Should show proper stellar distribution
- **Visible microlensing events** - Blink between peak and baseline images for F_peak > ~1.3
- **Proper scaling** - No white squares or completely black images

## 🆕 What's New (v2.0.0)

**For existing users:** All your workflows, parameter files, and scripts work exactly the same!

**New optional features:**
- **Input validation** - Catch configuration errors before running simulations
- **Comprehensive documentation** - Detailed guides at [Read the Docs](https://gulls.readthedocs.io)
- **Automated testing** - CI runs tests automatically
- **Better error messages** - Clear feedback when things go wrong
- **Version management** - Proper semantic versioning and releases

**See [CONTRIBUTING.md](CONTRIBUTING.md) for a gradual adoption guide.**

## What Gulls Does

Gulls simulates gravitational microlensing events with:
- **Single and binary sources/lenses**
- **Realistic observing conditions** (weather, seeing, detector effects)
- **Multiple observatories** and filter systems
- **Astrometric and photometric** signal generation
- **Detection statistics** and survey planning

## Executables

| Executable | Purpose |
|------------|---------|
| `gulls_std.x` | Standard microlensing simulation |
| `gulls_croin.x` | Crowding analysis and detection |
| `gullsFish.x` | Fisher matrix analysis |

## Input Files

Gulls requires several input files specified in your parameter file:

- **Sources**: Star catalogs with positions, magnitudes, and properties
- **Lenses**: Lens catalogs with masses, distances, and proper motions  
- **Observatories**: Telescope configurations, filters, and observing sequences
- **Starfields**: Background star distributions
- **Weather**: Observing conditions and weather profiles

See [PARAMETER_REFERENCE.md](PARAMETER_REFERENCE.md) for complete parameter documentation.

## Validating Your Input Files

Before running a simulation, validate your input catalogs to catch common issues:

```bash
python scripts/validate_inputs.py your_parameter_file.prm
```

This checks:
- Required column headers in source and lens catalogs
- Valid source/lens distance pairs (sources must be farther than lenses)
- Reasonable distance ranges (positive, < 50 kpc, not NaN/Inf)
- Binary source column requirements (when `MULTIPLE_SOURCES=1`)

The validation uses the same checks as CI, catching issues before you start long simulations.

## Examples

The `smoke_test/` directory contains working examples:
- `smoke_test/parameterfiles/` - Example parameter files
- `smoke_test/assets/` - Sample input catalogs and configurations

## Documentation

- [Parameter Reference](PARAMETER_REFERENCE.md) - Complete parameter documentation
- [Input File Formats](INPUT_FORMATS.md) - Catalog and configuration file specifications
- [Examples](examples/) - Working example configurations

## Contributing

1. Fork the repository
2. Create a feature branch
3. Run smoke tests: `python smoke_test/run_smoke_test.py`
4. Submit a pull request

## License

[License information]

## Citation

If you use Gulls in your research, please cite:
[Citation information]