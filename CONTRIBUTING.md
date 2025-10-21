# Contributing to Gulls

Thank you for your interest in contributing to Gulls! This document provides guidelines for contributing to the project.

## 🚀 Quick Start for Existing Developers

**Don't worry - everything you know still works!** The new features are **optional** and **additive**:

- ✅ **Your existing workflow** - Makefile, direct compilation, manual testing
- ✅ **Your existing tools** - Same executables, same parameter files  
- ✅ **Your existing scripts** - All your custom analysis code still works
- 🆕 **New optional features** - Validation, documentation, CI (use if you want)

### What's New (Optional)
- **Input validation** - `python3 scripts/validate_inputs.py your_file.prm` (catches errors early)
- **Better documentation** - See [Read the Docs](https://gulls.readthedocs.io) for detailed guides
- **Automated testing** - CI runs tests automatically on pull requests
- **Version management** - Automated version bumping and releases

### What Stays the Same
- **Same build process** - `make` still works exactly as before
- **Same executables** - `gulls_std.x`, `gulls_croin.x`, `gullsFish.x` unchanged
- **Same parameter files** - All your existing `.prm` files work identically
- **Same output format** - All output files are identical

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone <your-fork-url>`
3. Create a feature branch: `git checkout -b feature/your-feature-name`
4. Make your changes
5. Test your changes: `python3 smoke_test/run_smoke_test.py`
6. Submit a pull request

## Development Workflow

### Building Gulls

**Option 1: Traditional Makefile (unchanged)**
```bash
make clean
make
```

**Option 2: CMake (new, recommended)**
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

**Both methods produce identical executables!** Use whichever you prefer.

### Running Tests

**Option 1: Manual testing (unchanged)**
```bash
# Your existing workflow - test with your own parameter files
./bin/gulls_std.x -i your_file.prm -s 0 -f 0
```

**Option 2: Automated smoke tests (new, optional)**
```bash
# Run all smoke tests
python3 smoke_test/run_smoke_test.py

# Run CI subset (faster)
python3 smoke_test/run_smoke_test.py --ci

# Run specific test
python3 smoke_test/run_smoke_test.py --cases std-binary
```

**Option 3: Input validation (new, recommended)**
```bash
# Check your input files before running
python3 scripts/validate_inputs.py your_file.prm
```

## 🔄 Gradual Adoption Guide

**You don't need to change everything at once!** Here's how to gradually adopt new features:

### Phase 1: Just Try Validation (5 minutes)
```bash
# Before running your simulation, just try this:
python3 scripts/validate_inputs.py your_parameter_file.prm
```
This catches common errors early and saves debugging time.

If you run in to a new error, consider adding a check in the validation.

### Phase 2: Try CMake Build (10 minutes)
```bash
# Instead of 'make', try:
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
Same result, but with better dependency management.

### Phase 3: Use Documentation (when needed)
- Check [Read the Docs](https://gulls.readthedocs.io) when you have questions
- Look up parameter descriptions in the [Parameter Reference](https://gulls.readthedocs.io/en/latest/parameter_reference.html)

### Phase 4: Contribute Back (when ready)
- Use the new validation system for your contributions
- Follow the automated testing workflow
- Help improve documentation for others

**Remember: All old methods still work!** Adopt new features at your own pace.

### Code Style

- Follow existing C++ style conventions
- Use meaningful variable and function names
- Add comments for complex algorithms
- Update documentation for new features

## Types of Contributions

### Bug Reports

When reporting bugs, please include:
- GULLS version/commit
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Relevant log files

### Feature Requests

For new features, please:
- Check existing issues first
- Provide a clear description
- Explain the use case
- Consider backward compatibility

### Code Contributions

- Keep changes focused and atomic
- Add tests for new functionality
- Update documentation
- Ensure all tests pass

## Testing

All contributions must pass the smoke tests:

```bash
python3 smoke_test/run_smoke_test.py
```

The CI system will automatically run tests on pull requests.

## Documentation

When adding new features:
- Update relevant documentation files
- Update parameter reference if new parameters are added
- Add new smoke tests, where appropriate

## Pull Request Process

1. Ensure your branch is up to date with main/dev
2. Run smoke tests locally
3. Create a clear, descriptive pull request
4. Reference any related issues
5. Respond to review feedback promptly

## Questions?

Feel free to open an issue for questions about contributing or the codebase.
