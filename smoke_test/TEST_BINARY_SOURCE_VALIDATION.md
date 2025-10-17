# Binary Source Column Validation Tests

## Summary

Added validation to the smoke test that raises an error if `MULTIPLE_SOURCES=1` is set in the parameter file but the source catalog is missing required binary source columns.

## Answer to "Is 's' necessary?"

**No**, "s" is NOT a required input column. 

The variable `scomp_s` (projected separation in Einstein radii) is **calculated** from other columns in `buildEvent.cpp` line 474:
```cpp
Event->scomp_s.push_back(acomb/(Event->thE * Sources->data[sn][Sources->DIST]));
```

It's derived from:
- `combined_logP` (orbital period)
- `Mass` (stellar masses, via Kepler's 3rd law)
- `Dist` (distance)
- `thE` (Einstein radius)

## What Was Added

### 1. New Validation Function (`smoke_test/validation.py`)

Added `verify_binary_source_columns(params)` that:
- Checks if `MULTIPLE_SOURCES=1` in the parameter file
- If enabled, validates all source catalogs contain required columns
- Raises clear error messages listing missing columns

**Required Named Columns** (accessed via `datadict`):
- `Is_Binary` - binary status flag (0=single, 1=primary, 2+=companion)
- `ID` - unique identifier
- `primary_ID` - links companions to primaries
- `combined_logP` - log10 orbital period (days)

**Required Standard Columns** (accessed by position):
- `mul` - proper motion in l (mas/yr)
- `mub` - proper motion in b (mas/yr)
- `Mass` - stellar mass (solar masses)
- `Radius` - stellar radius (solar radii)
- `Dist` - distance (kpc)

### 2. Integration into Smoke Test Runner (`smoke_test/runner.py`)

The validation runs **before** executing any simulations, so it catches configuration errors early.

## How to Test

### Test 1: Valid Binary Source Catalog
```bash
cd /Users/malpas.1/Code/gulls_mp
python3 smoke_test/run_smoke_test.py --cases fish-binary
```
**Expected:** ✅ Test passes (catalog has all required columns)

### Test 2: Missing Column Detection

Create a test catalog missing `combined_logP`:

```bash
# Create test directory
mkdir -p smoke_test/assets/sources/test

# Copy the header but remove combined_logP
head -1 smoke_test/assets/sources/smoke_source_binary_catalog.dat | \
  sed 's/combined_logP//g' > smoke_test/assets/sources/test/broken_binary.dat

# Copy a few data rows
tail -n +2 smoke_test/assets/sources/smoke_source_binary_catalog.dat | \
  head -3 >> smoke_test/assets/sources/test/broken_binary.dat

# Create a sources list pointing to broken catalog
echo "0 0.0 0.0 0.0 broken_binary.dat" > smoke_test/assets/sources/test/broken.sources

# Create a test parameter file
cat smoke_test/parameterfiles/smoke_fish_binary.prm | \
  sed 's/SOURCE_LIST=.*/SOURCE_LIST=test\/broken.sources/' > \
  smoke_test/parameterfiles/test_broken_binary.prm
```

Then run:
```bash
python3 smoke_test/run_smoke_test.py --cases test-broken-binary
```

**Expected:** ❌ Test fails with clear error:
```
Binary source validation failed for test-broken-binary:
 - MULTIPLE_SOURCES=1 but source catalog broken_binary.dat 
   is missing required binary source columns: combined_logP
```

## Example Error Messages

### Missing Binary-Specific Columns
```
MULTIPLE_SOURCES=1 but source catalog smoke_source_binary_catalog.dat 
is missing required binary source columns: Is_Binary, combined_logP, primary_ID
```

### Missing Standard Columns
```
MULTIPLE_SOURCES=1 but source catalog smoke_source_binary_catalog.dat 
is missing required standard columns: Dist, Mass, Radius
```

### Missing SOURCE_DIR or SOURCE_LIST
```
MULTIPLE_SOURCES=1 requires SOURCE_DIR and SOURCE_LIST to be set
```

## Code References

### Where Columns Are Used

**In `buildEvent.cpp`:**
- Line 388: `Is_Binary` - identifies binary systems
- Line 398: `primary_ID` & `ID` - links companions to primaries
- Line 467: `RADIUS`, `DIST` - calculates source size in Einstein radii
- Line 468: `combined_logP` - orbital period
- Line 469-470: `Mass` - stellar masses for Kepler's law
- Line 474: Calculates `scomp_s` from period, masses, distance

**In `pllxLightcurveGenerator.cpp`:**
- Lines 107-112: Uses `scomp_s`, `scomp_phase`, `scomp_inc`, `scomp_alpha` (all calculated from the required columns)
- Line 112: Uses `scomp_rs` (calculated from `RADIUS` and `DIST`)
- Line 117: Uses `scomp_fsofs1` (calculated from magnitudes)

## Benefits

1. **Early Error Detection**: Catches missing columns before running expensive simulations
2. **Clear Error Messages**: Tells users exactly which columns are missing
3. **Prevents Silent Failures**: No more mysterious crashes due to missing binary source data
4. **Documentation**: The validation itself serves as documentation of required columns

