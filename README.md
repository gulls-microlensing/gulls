# Gulls: Gravitational Microlensing Simulation Suite

Gulls is a comprehensive simulation framework for modeling gravitational microlensing events with realistic observing conditions.

## Quick Start

### Prerequisites
- **CMake** 3.10 or later
- **C++17** compatible compiler (GCC, Clang, or MSVC)
- **Fortran** compiler (gfortran recommended)
- **Python 3.7+** (for smoke tests and validation)

### Build
```bash
git clone <repository-url>
cd gulls_push
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Validate Your Input Catalogs (Recommended)
```bash
python3 scripts/validate_inputs.py parameter_file.prm
```

### Run a Simulation
```bash
./build/bin/gulls_std.x parameter_file.prm
```

### Test the Installation
```bash
python3 smoke_test/run_smoke_test.py --ci
```

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
python3 scripts/validate_inputs.py your_parameter_file.prm
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
3. Run smoke tests: `python3 smoke_test/run_smoke_test.py`
4. Submit a pull request

## License

[License information]

## Citation

If you use Gulls in your research, please cite:
[Citation information]