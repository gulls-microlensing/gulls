# Gulls v2.0.0 Release Notes

**Release Date:** October 2025

## Major Release - Documentation and Testing Overhaul

This is a major release that significantly improves the usability and maintainability of Gulls for the broader microlensing community.

## What's New

### 📚 Comprehensive Documentation
- **Complete documentation system** with Sphinx/Read the Docs integration
- **Input file format specifications** for all file types (catalogs, observatories, sequences)
- **Parameter reference** with detailed descriptions of all configuration options
- **Installation guides** with both CMake and traditional Makefile approaches
- **Troubleshooting guides** for common issues

### 🔍 Input Validation System
- **Comprehensive validation** of input catalogs and configuration files
- **File existence checks** - ensures all referenced files exist
- **Data format validation** - checks column headers, data types, ranges
- **Compatibility checks** - validates source/lens distance relationships
- **Binary source validation** - ensures required columns when using multiple sources
- **Clear error messages** with actionable fixes

### 🧪 Testing and CI/CD
- **Smoke test suite** with automated testing of core functionality
- **GitHub Actions CI** with multi-platform testing (Ubuntu, macOS)
- **Automated validation** in CI pipeline
- **Example configurations** for different simulation types
- **Visual validation outputs** with smoke test output figures in Release Notes

### 🛠️ Developer Experience
- **CMake build system** alongside traditional Makefile
- **Contributing guidelines** for community contributions
- **Version management** with automated bumping and release workflows
- **Code quality improvements** with better error handling
- **Simplified release workflow**: patch > edit changelog > release

### 🎯 PSF Generation Improvements
- **Smoke tests generate PSF files on-demand** instead of requiring pre-committed files

## Breaking Changes

- **Buffer size fixes** - Fixed potential buffer overflows in path handling
- **Input validation** - Stricter validation may catch previously ignored configuration errors
- **Documentation structure** - New documentation format (RST instead of markdown)

## Migration Guide

### For Existing Users
1. **Update your build process** - CMake is now recommended over Makefile
2. **Validate your input files** - Run `python scripts/validate_inputs.py your_file.prm` before simulations
3. **Check documentation** - New comprehensive guides available at [Read the Docs](https://gulls.readthedocs.io)

### For Developers
1. **Use the new validation system** - Add validation for new error conditions
2. **Follow contributing guidelines** - See `CONTRIBUTING.md` for development workflow
3. **Update version numbers** - Use `python scripts/bump_version.py` for releases

## Technical Improvements

### Bug Fixes
- Fixed infinite loop in random number generation CI stub
- Fixed uninitialized memory issues in binary source calculations  
- Fixed off-by-one errors in catalog parsing
- Fixed buffer overflows in file path construction (required change for successful CI runs)
- Fixed failure to build docs in the release workflow on GitHub
- Fixed PSF generation in CI environments (removed hardcoded local machine paths)
- Fixed PSF file size issues (now generates proper 68MB files with subpixel sampling)
- Fixed simulation crashes due to missing or malformed PSF files

### Performance
- Improved error handling and user feedback
- Optimized validation routines
- Smart PSF caching reduces redundant file generation

### Security
- Fixed potential buffer overflows
- Improved input sanitization
- Better error handling to prevent crashes

## Community Impact

This release makes Gulls significantly more accessible to the broader microlensing community:

- **Easier installation** with better dependency management
- **Clear documentation** for new users
- **Robust validation** prevents common configuration errors
- **Professional development workflow** for contributors
- **Automated testing** ensures reliability
- **Improved CI/CD** with proper documentation builds and PSF handling

## Acknowledgments

This release represents a major community effort to improve Gulls' usability and maintainability. Special thanks to all contributors who helped with documentation, testing, and code improvements.

## Getting Started

1. **Install Gulls** - See the [Installation Guide](https://gulls.readthedocs.io/en/latest/install_gulls.html)
2. **Validate your inputs** - Use `python scripts/validate_inputs.py your_file.prm`
3. **Run simulations** - See the [Running Guide](https://gulls.readthedocs.io/en/latest/run_simulations.html)
4. **Get help** - Check the [Troubleshooting Guide](https://gulls.readthedocs.io/en/latest/basic_troubleshooting.html)

## What's Included

- **Source code**: Complete Gulls source with CMake build system
- **Binaries**: Linux executables (GSL fallbacks - testing only)
- **Documentation**: Built HTML documentation
- **Smoke test plots**: Visual proof that the release works

## Full Changelog

See [CHANGELOG.md](CHANGELOG.md) for the complete list of changes.

---

**Previous Release:** v1.0.0 (2013-2025)