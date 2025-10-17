# Debug vs Release Build Comparison: Binary Source Bug

**Date:** October 15, 2025  
**Test:** Binary source simulation with `smoke_binary.sources` catalog  
**Seed:** 12345  
**LC_TIMEOUT:** 600 seconds

## Executive Summary

**BOTH Debug and Release builds exhibit the bug** (`inf` values for `Source2_s` and `Source2_rho`), but with different failure modes:

- **Debug build:** VBM times out at 600s limit
- **Release build:** VBM completes in 365s but produces invalid output (`nan` flux values)

## Detailed Results

### Debug Build (-g)

```
Execution time: 10m 2s (602 seconds)
Total: 602.727 sec
LightcurveGen: 601.986 sec
Result: "Discarding event 0 (Failed lightcurve generation)"
Output: No .lc files generated
Source2_s: (not written to output)
Source2_rho: (not written to output)
```

**Behavior:** Hit the `LC_TIMEOUT=600.0` limit and was discarded.

### Release Build (-O3)

```
Execution time: 6m 5s (365 seconds)
Total: 365.379 sec  
LightcurveGen: 364.825 sec
Result: "All smoke test cases completed successfully"
Output: Generated 1661-line .lc file and plots
Source2_s: inf
Source2_rho: inf
Lightcurve flux: nan (invalid)
```

**Behavior:** VBM completed faster (likely due to compiler optimizations), but produced invalid results.

Sample lightcurve output:
```
t=0.00289351851852 mag=nan
t=0.124421296296 mag=nan
t=0.245949074074 mag=nan
```

## Key Findings

1. **The bug exists in both builds** - It's not a debug-specific initialization issue
2. **Compiler optimizations affect runtime but not correctness** - Release mode runs ~40% faster but still produces `inf` values
3. **`inf` companion parameters lead to invalid lightcurves** - VBM doesn't handle `inf` gracefully
4. **Matt's "successful" runs likely saw the Release behavior** - Quick completion with plots, but invalid flux values that may not have been checked

## Root Cause Confirmation

The bug is definitively in `buildEvent.cpp` lines 472-480:
```cpp
Event->scomp_s = acomb/(Event->thE*dist);      // Line 473 - thE uninitialized!
Event->scomp_rs = rscomb/Event->thE;           // Line 480 - thE uninitialized!
```

Where `Event->thE` is only assigned later at line 516:
```cpp
Event->thE = Event->rE/dist;
```

With uninitialized `thE`, division produces `inf`, regardless of build type.

## Implications for Matt's Tests

Matt's tests likely succeeded in the sense that they:
- ✅ Ran to completion (Release build)
- ✅ Generated plots
- ❌ But produced invalid flux measurements (`nan` values)

If Matt didn't check the actual flux values in the lightcurves, the runs would have appeared successful.

## Recommendation

**Apply the fix** from `gulls_mp` branch to move the `thE` calculation before its use:

```cpp
// Calculate thE first (around line 466, BEFORE companion calculations)
Event->thE = Event->rE/dist;

// Then use it for companion parameters
if(Source->scomp>0){
    Event->scomp_s = acomb/(Event->thE*dist);
    Event->scomp_rs = rscomb/Event->thE;
    // ... rest of companion logic
}
```

This is a clear logic error, not a configuration or catalog issue.

