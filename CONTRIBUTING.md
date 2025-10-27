# Contributing to Gulls

Thank you for your interest in contributing to Gulls! This document provides guidelines for contributing to the project.

## Quick Start for Existing Developers

**Don't worry - everything you know still works!** The new features are **optional** and **additive**:

- ✅ **Your existing workflow** - Makefile, direct compilation, manual testing
- ✅ **Your existing tools** - Same executables, same parameter files  
- ✅ **Your existing scripts** - All your custom analysis code still works
- 🆕 **New optional features** - Validation, documentation, CI (use if you want)

### What's New (Optional)
- **Input validation** - `python scripts/validate_inputs.py your_file.prm` (catches errors early)
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
5. Test your changes: `python smoke_test/run_smoke_test.py`
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
python smoke_test/run_smoke_test.py

# smoke_test help text (lists options and cases)
python smoke_test/run_smoke_test.py --help

# Run CI subset (faster)
python smoke_test/run_smoke_test.py --ci

# Run specific test
python smoke_test/run_smoke_test.py --cases std-binary
```

**Option 3: Input validation (new, recommended)**
```bash
# Check your input files before running
python scripts/validate_inputs.py your_file.prm
```

## 🔄 Gradual Adoption Guide

**You don't need to change everything at once!** Here's how to gradually adopt new features:

### Phase 1: Just Try Validation (5 minutes)
```bash
# Before running your simulation, just try this:
python scripts/validate_inputs.py your_parameter_file.prm
```
This catches common errors early and saves debugging time.

If you run in to a new error, consider adding a check in the validation.

### Phase 2: Try CMake Build (10 minutes)
```bash
# Instead of 'make', try:
cmake -S . -B build
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

## 🔍 Adding New Validation Tests

One of the best ways to contribute is by adding validation for new error conditions you encounter. This helps prevent others from hitting the same issues.

### When to Add Validation

Add validation when you encounter:
- ✅ **Configuration errors** that cause crashes or hangs
- ✅ **Data format issues** that lead to incorrect results  
- ✅ **Missing files** that cause cryptic error messages
- ✅ **Invalid parameter combinations** that cause undefined behavior
- ❌ **Performance issues** (validation is for correctness, not speed)

### How to Add Validation

**Step 1: Add validation function to `smoke_test/validation.py`**

```python
def verify_your_new_condition(params: Dict[str, str]) -> None:
    """Check that [your condition] is satisfied."""
    # Load relevant data
    some_file = _resolve_param_path(params.get("SOME_FILE"), "some file")
    
    if not some_file.exists():
        raise SmokeTestError(f"Some file not found: {some_file}")
    
    # Check your condition
    if condition_violated:
        raise SmokeTestError(
            "Clear error message explaining what's wrong "
            "and how to fix it"
        )
```

**Step 2: Export it in `__all__` at the bottom of `validation.py`**

```python
__all__ = [
    # ... existing functions ...
    "verify_your_new_condition",
]
```

**Step 3: Call it in `smoke_test/runner.py`**

Add to the validation loop around line 154:
```python
verify_your_new_condition(case.params)
```

**Step 4: Test your validation**

```bash
# Test with a bad file to make sure it catches the error
python scripts/validate_inputs.py your_test_file.prm
```

### Example: Adding Weather File Validation

```python
def verify_weather_coverage(params: Dict[str, str]) -> None:
    """Verify weather file covers the full simulation duration."""
    sim_length = float(params.get("SIMULATION_LENGTH", "0"))
    weather_dir = _resolve_param_path(params.get("WEATHER_PROFILE_DIR"), "weather directory")
    
    # Check weather files exist and cover simulation length
    weather_files = list(weather_dir.glob("*.weather"))
    if not weather_files:
        raise SmokeTestError(f"No weather files found in {weather_dir}")
    
    # Parse weather file to check max time
    max_time = 0.0
    for weather_file in weather_files:
        # ... parse file and check coverage ...
    
    if max_time < sim_length:
        raise SmokeTestError(
            f"Weather file only covers {max_time:.2f} days "
            f"but SIMULATION_LENGTH is {sim_length:.2f} days.\n"
            f"Use scripts/make_weather.py to generate longer weather profiles."
        )
```

### Benefits of Adding Validation

- **Prevents others** from hitting the same error
- **Clear error messages** instead of cryptic crashes
- **Documentation** of what can go wrong
- **Automated testing** catches regressions
- **Community contribution** that helps everyone

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
python smoke_test/run_smoke_test.py
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
