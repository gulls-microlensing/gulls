# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-10-XX

### Added
- Comprehensive documentation system with Sphinx/Read the Docs
- Input file format specifications and parameter reference
- Validation system for input catalogs and configuration files
- CI/CD pipeline with GitHub Actions
- Smoke test suite for automated testing
- CMake build system alongside traditional Makefile
- Contributing guidelines and development workflow
- Example configurations and troubleshooting guides

### Changed
- Improved error handling and user feedback
- Enhanced input validation with detailed error messages
- Streamlined installation process with better dependency management
- Updated documentation from legacy format to modern RST

### Fixed
- Buffer overflow issues in path handling
- Infinite loop bugs in random number generation
- Uninitialized memory issues in binary source calculations
- Off-by-one errors in catalog parsing

### Security
- Fixed potential buffer overflows in file path construction
- Improved input validation to prevent malformed data crashes

## [1.0.0] - 2013-2025

### Added
- Core microlensing simulation framework
- Support for single and binary sources/lenses
- Realistic observing conditions and detector effects
- Multiple observatory and filter system support
- Photometric signal generation
- Lens orbital motion
- Support for low-level VBMicrolensing functions
- Detection statistics and survey planning tools
- Many other features, executables, and documentation evolutions.

### Original Development
- Initial implementation by Matthew Penny and collaborators
- Published in Penny et al. (2013, 2014, 2019)
- Core algorithms for gravitational microlensing event simulation
