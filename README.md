# gulls

A microlensing simulator optimized for space-based microlensing
surveys, but also supporting ground-based observatory simulations.   

## Requirements

1. A C++ compiler with C++17 support (e.g. `g++`, `clang++`, `icpx`)
1. A Fortran compiler (e.g. `gfortran`, `ifx`)
1. [CMake ≥ 3.20](https://cmake.org/download/)
1. [GNU Scientific Library](https://www.gnu.org/software/gsl/)
1. [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
1. `ESPL.tbl` from the
   [VBMicrolensing data directory](https://github.com/valboz/VBMicrolensing/tree/master/VBMicrolensing/data)

> **Heads-up:** the repository now vendors the latest
> `VBMicrolensingLibrary.cpp`/`.h`, so you no longer need a pre-built
> `libVBB.a` to link against. You only need the `ESPL.tbl` lookup table.

## Building with CMake (recommended)

1. Install the requirements above. On macOS you can use Homebrew:
   ```bash
   brew install cmake gsl cfitsio gcc
   ```
   On Linux, your package manager typically provides the same
   dependencies.
1. Clone this repository (and optionally place the
   [`VBMicrolensing`](https://github.com/valboz/VBMicrolensing)
   repository alongside it so you can copy `ESPL.tbl`).
1. Copy `ESPL.tbl` into `gulls_mp/src/`.
1. Configure and build:
   ```bash
   cmake -S gulls_mp -B gulls_mp/build
   cmake --build gulls_mp/build
   ```
   This produces `gulls_std`, `gulls_croin`, and `gullsFish` in
   `gulls_mp/build/bin/`.
1. (Optional) Install the binaries anywhere you like with
   `cmake --install build --prefix <path>`.

## Running the executables

After building, the executables are located in `build/bin/`:

```bash
# Run directly from the project root directory
./build/bin/gulls_std <parameter_file> [options]
./build/bin/gulls_croin <parameter_file> [options]
./build/bin/gullsFish <parameter_file> [options]
```

Or add the build directory to your PATH for easier access:

```bash
export PATH="$PWD/build/bin:$PATH"
gulls_std <parameter_file> [options]
```

Use the `-d` flag for debug output (repeat for more verbosity: `-d`, `-dd`, `-ddd`).

### Selecting a build type

By default the project configures in `Release` mode. To switch to
`Debug` (with symbols and runtime checks) pass

```bash
cmake -S gulls_mp -B build -DCMAKE_BUILD_TYPE=Debug
```

### Compiler warnings

By default, compiler warnings are **disabled** for a cleaner build output. 
To enable warnings (useful when fixing code issues):

```bash
cmake -B build -DENABLE_WARNINGS=ON && cmake --build build
```

To disable warnings again:

```bash
cmake -B build -DENABLE_WARNINGS=OFF && cmake --build build
```

For quick rebuilds with the current warning setting, just use:

```bash
cmake --build build
```

### Legacy Makefile workflow

The historical Makefiles remain in `src/` for anyone who relies on the
old process. They still expect a lot of manual setup (absolute paths,
Intel compilers, etc.) and are **no longer recommended**. If you must
use them, run `./configure.sh` to rewrite hard-coded paths and then
`make <target>` inside `src/`.

> Known gaps: a handful of older targets reference Fortran sources such
> as `findTrack.f`, `magTrack.f`, `extend.f`, and `readData.f` that are
> not present in the repository. Only the `gulls_std`, `gulls_croin`,
> and `gullsFish` executables currently build successfully.

## (Incomplete) Checklist/Troubleshooting for regular gulls runs

- When preparing a run it's a good idea to do a short run with 
    ```
    PRETTY_PICS=1
    PRETTY_PICS_DIMENSIONS=128 #Or bigger if you'd like
    OUTPUT_LC=1 #If this is less than 1, it represents a probability of output
    OUTPUT_IMAGES=1
    OUTPUT_ONALL=1
    ```
    This should produce the maximum amount of diagnostic
    output. `OUTPUT_LC` determines the frequency with which lightcurves
    are output; `OUTPUT_ONALL` will output every event generated, but if
    set to zero, `OUTPUT_ONDET` will only output when a detection occurs
    (definition of detection depends on the type of run). `OUTPUT_IMAGES`
    saves images each time a lightcurve is output - one at peak
    magnification, one at baseline for each filter. Inspect both the
    images and the lightcurves. .det.lc files should show some level of
    detectable microlensing event in the lightcurve, .all.lc may
    not. Images should show a starfield, and depending on the peak
    magnification, blinking between peak and base images should show the
    microlensing event (requires Fpeak>~1.3 to be easily visible).

    - If the image is a white square, first make sure you select zscale in
ds9, and if that doesn't reveal stars it might indicate the presence
of bad magnitudes in the star, source, or lens catalog(s), or it's
possible that you just got unlucky and had a super-bright star land in
the image - the larger the PRETTY_PIC, the less likely it is that such
a star will blow out the image. Try generating a very large pretty pic
and see how many super-bright stars there are (they will look like
white squares) -- too many of these and its probably worth looking at
the color magnitude diagrams of the star lists to see that there's
nothing weird in there. 
- You can run with different levels of debug information by repeating -d flags
    - If you get an message about no valid stars, check that you have
	`NFILTERS` set correctly. 

