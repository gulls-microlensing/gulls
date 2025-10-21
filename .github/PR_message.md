# Add CI/CD, Documentation, and Testing Infrastructure

## Summary
This PR adds modern development infrastructure to Gulls while maintaining full backward compatibility. All existing workflows, parameter files, and scripts continue to work unchanged.

## What's New

### 🔧 **CI/CD Pipeline**
- **GitHub Actions CI** - Automated testing on Ubuntu/macOS
- **Release workflow** - Automated releases when you push version tags
- **Documentation builds** - Sphinx/Read the Docs integration

### 📚 **Documentation System**
- **Comprehensive guides** - Installation, input formats, parameter reference
- **Contributing guidelines** - Clear path for community contributions
- **Troubleshooting guides** - Common issues and solutions

### 🧪 **Testing & Validation**
- **Smoke test suite** - Automated testing of core functionality
- **Input validation** - `python scripts/validate_inputs.py your_file.prm`
- **Version management** - `python scripts/bump_version.py patch/minor/major`

### 🔄 **Version Management**
- **Semantic versioning** - Proper v2.0.0 instead of hardcoded dates
- **Automated bumping** - Updates version in both `gulls.cpp` and documentation
- **Release automation** - Creates GitHub releases with archives

## Usage Instructions

### For Existing Users
**Nothing changes!** Your existing workflow, parameter files, and scripts work identically.

### For New Features
```bash
# Validate inputs before running
python scripts/validate_inputs.py your_file.prm

# Bump version and create release
python scripts/bump_version.py patch    # 2.0.0 -> 2.0.1
# Edit CHANGELOG.md with your changes
python scripts/bump_version.py release  # Commit, tag, and push automatically
```

### For Matt's Workflow
```bash
# 1. Make your changes
# 2. Bump version
python scripts/bump_version.py patch

# 3. Edit CHANGELOG.md (add your changes)
# 4. Create release (handles everything automatically)
python scripts/bump_version.py release
```

**That's it!** The script handles:
- ✅ Commits changes (with smart unstaged change detection)
- ✅ Creates and pushes tags
- ✅ Triggers release workflow automatically
- ✅ Handles existing tags gracefully

## Workflow Triggers

### CI Workflow (`.github/workflows/test.yml`)
- **Triggers**: Push to any branch, pull requests
- **What it does**: Builds Gulls, runs smoke tests, validates inputs
- **Manual trigger**: Go to Actions tab → "Test" → "Run workflow"

### Release Workflow (`.github/workflows/release.yml`)
- **Triggers**: Push tags matching `v*` (e.g., `v2.0.0`, `v2.0.1`)
- **What it does**: Builds, tests, creates GitHub release with source/binary archives
- **How to trigger**: 
  ```bash
  # Manual (old way)
  git tag v2.0.0
  git push origin v2.0.0
  
  # Automated (new way)
  python scripts/bump_version.py release
  ```

### Documentation Workflow (`.github/workflows/docs.yml`)
- **Triggers**: Push to `main` branch
- **What it does**: Builds Sphinx documentation
- **Note**: May be redundant with separate `gulls-microlensing.github.io` repo

## Open Questions for Review

1. **Documentation hosting**: Currently builds docs in this repo, but `gulls-microlensing.github.io` exists separately. Should we:
   - Make `gulls-microlensing.github.io` a submodule of this repo?
   - Remove the docs workflow from this repo?
   - Keep both (redundant but safe)?

2. **Version strategy**: The version bumping script creates Git tags. Do you want to:
   - Use it for official releases?
   - Keep manual version management?

## Technical Changes
- **Buffer size fixes** - Required for CI environment (long paths)
- **Version updates** - v2.0.0 with proper date (October 2025)
- **Stub implementations** - GSL fallbacks for CI (Numerical Recipes preferred for production)

## Files Changed
- Added: CI workflows, documentation, validation scripts, version management
- Removed: Build artifacts (hundreds of files - makes diff look larger than it is)
- Modified: Buffer sizes, version numbers, added stubs