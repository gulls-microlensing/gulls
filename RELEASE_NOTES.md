# Gulls v2.0.1 Release Notes

**Release Date:** 2025-10-22

## Major Release


### Added
- Missing docs build requirements in the `environment.yml`
- Smoke test output figures in Release Notes
- PSF generation utilities (`generateMoffat`, `txt2fits_custom`, `precompute_psf`) to CMake build
- On-demand PSF file generation for smoke tests (no more 68MB files in git)
- Smart PSF caching system that reuses existing files when available

### Changed
- Simplified release workflow: patch > edit changelog > release
- PSF generation now uses proper subpixel sampling (9×9 = 81 variations)
- Smoke tests generate PSF files on-demand instead of requiring pre-committed files

### Fixed
- Failure to build docs in the release workflow on GitHub
- PSF generation in CI environments (removed hardcoded local machine paths)
- PSF file size issues (now generates proper 68MB files with subpixel sampling)
- Simulation crashes due to missing or malformed PSF files
## What's New

This release includes the following changes:

## What's Included

- **Source code**: Complete Gulls source with CMake build system
- **Binaries**: Linux executables (GSL fallbacks - testing only)
- **Documentation**: Built HTML documentation
- **Smoke test plots**: Visual proof that the release works

## Getting Started

1. **Install Gulls** - See the [Installation Guide](https://gulls.readthedocs.io/en/latest/install_gulls.html)
2. **Validate your inputs** - Use `python scripts/validate_inputs.py your_file.prm`
3. **Run simulations** - See the [Running Guide](https://gulls.readthedocs.io/en/latest/run_simulations.html)

## Full Changelog

See [CHANGELOG.md](CHANGELOG.md) for the complete list of changes.

---

**Previous Release:** v2.0.0