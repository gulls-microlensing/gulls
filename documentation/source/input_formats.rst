Input File Formats
==================

This document describes the required formats for Gulls input files.

.. tip::
   Validate your input catalogs before running simulations:
   
   .. code-block:: bash
   
      python3 scripts/validate_inputs.py your_parameter_file.prm
   
   This checks for required columns, valid source/lens pairs, reasonable
   distance ranges, and binary source requirements.

Source Catalogs
---------------

Source catalogs contain information about the stars that can be lensed.

Catalog Format
~~~~~~~~~~~~~~

Source catalogs are whitespace-delimited text files with:

1. **Header line**: Column names (must match exactly as shown below)
2. **Data lines**: One star per line

.. important::
   The **first N columns** must be **photometric magnitudes** in various filters, where N = ``NFILTERS``
   parameter in your ``.prm`` file. These filter names should correspond to filters used in your
   observatory sequence files.

**Example filter columns** (first 32 columns for NFILTERS=32):

- Roman filters: ``R062``, ``Z087``, ``Y106``, ``J129``, ``W146``, ``H158``, ``F184``, ``K213``
- 2MASS: ``2MASS_Ks``, ``2MASS_J``, ``2MASS_H``
- Bessell: ``Bessell_U``, ``Bessell_B``, ``Bessell_V``, ``Bessell_I``, ``Bessell_R``
- Other systems: ``Kepler_Kp``, ``TESS``, ``DECam_*``, ``Gaia_*``, ``VISTA_*``

Required Columns (after magnitude columns)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| mul            | float  | Proper motion in l             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| mub            | float  | Proper motion in b             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| Mass           | float  | Stellar mass                   | M☉       |
+----------------+--------+--------------------------------+----------+
| Radius         | float  | Stellar radius                 | R☉       |
+----------------+--------+--------------------------------+----------+
| Dist           | float  | Distance to source             | kpc      |
+----------------+--------+--------------------------------+----------+

Optional Columns
~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| RA2000.0       | float  | Right ascension (J2000)        | degrees  |
+----------------+--------+--------------------------------+----------+
| DEC2000.0      | float  | Declination (J2000)            | degrees  |
+----------------+--------+--------------------------------+----------+
| Vr             | float  | Radial velocity                | km/s     |
+----------------+--------+--------------------------------+----------+
| [Fe/H]         | float  | Metallicity                    | dex      |
+----------------+--------+--------------------------------+----------+
| Teff           | float  | Effective temperature          | K        |
+----------------+--------+--------------------------------+----------+
| logg           | float  | Surface gravity                | cgs      |
+----------------+--------+--------------------------------+----------+

Binary Source Columns
~~~~~~~~~~~~~~~~~~~~~

When ``MULTIPLE_SOURCES=1``, additional columns are required:

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| Is_Binary      | int    | Binary flag (0=single,         | -        |
|                |        | 1=primary, 2=companion)        |          |
+----------------+--------+--------------------------------+----------+
| ID             | int    | Unique source identifier       | -        |
+----------------+--------+--------------------------------+----------+
| primary_ID     | int    | ID of primary star             | -        |
|                |        | (for companions)               |          |
+----------------+--------+--------------------------------+----------+
| combined_logP  | float  | Log10 of orbital period        | days     |
+----------------+--------+--------------------------------+----------+

Source List File (``.sources``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your parameter file specifies ``SOURCE_DIR`` and ``SOURCE_LIST``. The ``SOURCE_LIST`` file
(e.g., ``my_sources.sources``) maps fields to catalog files:

**Format**: ``field_number l b l_width b_width catalog_filename``

+----------------+--------+------------------------------------------+
| Column         | Type   | Description                              |
+================+========+==========================================+
| field_number   | int    | Field identifier (0, 1, 2, ...)          |
+----------------+--------+------------------------------------------+
| l              | float  | Galactic longitude of field center (°)   |
+----------------+--------+------------------------------------------+
| b              | float  | Galactic latitude of field center (°)    |
+----------------+--------+------------------------------------------+
| l_width        | float  | Field width in l direction (degrees)     |
+----------------+--------+------------------------------------------+
| b_width        | float  | Field height in b direction (degrees)    |
+----------------+--------+------------------------------------------+
| catalog_file   | string | Name of catalog file (in SOURCE_DIR)     |
+----------------+--------+------------------------------------------+

**Example** (``smoke.sources``):

.. code-block:: text

   0 0.0 0.0 0.01 0.01 smoke_source_catalog.dat

This defines field 0 at Galactic coordinates (l, b) = (0.0°, 0.0°) with a 0.01° × 0.01° size,
using the catalog ``SOURCE_DIR/smoke_source_catalog.dat``.

Lens Catalogs
-------------

Lens catalogs have the same format as source catalogs: magnitude columns first, then required physical properties.

Required Columns (after magnitude columns)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| mul            | float  | Proper motion in l             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| mub            | float  | Proper motion in b             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| Mass           | float  | Lens mass                      | M☉       |
+----------------+--------+--------------------------------+----------+
| Dist           | float  | Distance to lens               | kpc      |
+----------------+--------+--------------------------------+----------+

.. note::
   Lens catalogs typically have no magnitude limit, as mass (not brightness) is the
   primary factor for lensing. The number of lenses needed depends on field pointing
   and solid angle (e.g., ~10⁴ for a crowded Galactic bulge field).

Lens List File (``.lenses``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Same format as source list files:

**Format**: ``field_number l b l_width b_width catalog_filename``

**Example** (``smoke.lenses``):

.. code-block:: text

   0 0.0 0.0 0.01 0.01 smoke_lens_catalog.dat

Starfield Catalogs
------------------

Starfield catalogs contain background stars (non-lensing, non-source stars) for realistic image generation.

Catalog Format
~~~~~~~~~~~~~~

Starfield catalogs have the same format as source/lens catalogs: magnitude columns first, then physical properties.
The catalog size and area/weight parameters (in the ``.starfields`` list file) should be tuned iteratively
to achieve the desired number of stars sampled per image for your field and detector configuration.

**Magnitude categories** (for efficient sampling):

- Bright: H ≤ 15 mag
- Moderate 1: 15 < H ≤ 20 mag  
- Moderate 2: 20 < H ≤ 25 mag
- Faint: H > 25 mag

Starfield List File (``.starfields``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Different format** from source/lens lists - includes an area/weight column for sampling:

**Format**: ``field_number l b area_or_weight catalog_filename``

+------------------+--------+------------------------------------------------------+
| Column           | Type   | Description                                          |
+==================+========+======================================================+
| field_number     | int    | Field identifier (0, 1, 2, ...)                      |
+------------------+--------+------------------------------------------------------+
| l                | float  | Galactic longitude of field center (°)               |
+------------------+--------+------------------------------------------------------+
| b                | float  | Galactic latitude of field center (°)                |
+------------------+--------+------------------------------------------------------+
| area_or_weight   | float  | Sampling weight (for brightness level stratification)|
+------------------+--------+------------------------------------------------------+
| catalog_file     | string | Name of starfield catalog file (in STARFIELD_DIR)    |
+------------------+--------+------------------------------------------------------+

**Example** (``smoke.starfields``):

.. code-block:: text

   0 0 1.0 smoke_field0_level0.starfield

.. note::
   The area/weight column enables proper sampling for very bright stars, which are rare
   but important for realistic images. Different magnitude levels may have different weights.

Observatory Configuration
-------------------------

Observatory files define telescope and detector properties.

Required Files
~~~~~~~~~~~~~~

- ``observatory.observatory`` - Main observatory configuration
- ``observatory.list`` - List of observatory files
- ``observatory.sequence`` - Observing sequence
- ``observatory.throughput`` - Filter throughput curves
- ``observatory.detector`` - Detector properties
- ``observatory.centres`` - Observatory locations

Observing Sequence Files (``.sequence``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequence files define the observing cadence and time allocation for each field.

**Format**:

.. code-block:: text

   #Key:
   #Nstack +ve, Texp +ve  (Image being taken by this instrument)
   #Nstack -ve, Texp +ve  (Image being taken by other instrument)
   #Nstack -ve, Texp -ve  (Overhead: slew, readout, filter change, etc.)

   #Field  Nstack  Texp    Sum     Description
   BEGIN_REPEAT 15
   0       -1      66.88   66.88   field 0 - F146 exposure (other observatory)
   0       -1     -661.628 728.51  Other fields (overhead)
   END_REPEAT

   0        1      66.88   66.88   field 0 - F087 exposure (THIS observatory)
   0       -1     -661.628 728.51  Other fields (overhead)

**Column meanings**:

+----------+--------+--------------------------------------------------------------+
| Column   | Type   | Description                                                  |
+==========+========+==============================================================+
| Field    | int    | Field number being observed                                  |
+----------+--------+--------------------------------------------------------------+
| Nstack   | int    | +: this observatory observes; -: time passing/overhead       |
+----------+--------+--------------------------------------------------------------+
| Texp     | float  | +: exposure time (s); -: overhead time (s)                   |
+----------+--------+--------------------------------------------------------------+
| Sum      | float  | Running total time (for reference/cadence calculation)       |
+----------+--------+--------------------------------------------------------------+
| Descrip. | string | Comment describing this sequence step                        |
+----------+--------+--------------------------------------------------------------+

**Repeat blocks**: Use ``BEGIN_REPEAT N`` / ``END_REPEAT`` to loop sequence sections N times.

.. important::
   **Simulation cadence** = sum of all sequence times (including repeats).
   
   Additional readout time between exposures in a stack is specified by ``READ_OHEAD``
   in the ``.observatory`` file, not in the sequence.

**Example interpretation** (from smoke.sequence):

.. code-block:: text

   0   1   500.0   500.0   field 0 - F184 daily exposure
   0  -1  10000.0 10000.0  other instrument

This means:
- Field 0: 500s exposure by this observatory (Nstack=1, Texp=500)
- Then 10,000s passes while other instruments observe (Nstack=-1, Texp=10000)
- **Total cadence**: 10,500 seconds between field 0 observations

**Complex example with repeats**:

From ``Roman_F146_overguide_6hcc.sequence``:

- 15 × (F146 exposure by others + overhead) = ~3h
- 1 × (F087 exposure + overhead) = ~12 min
- 15 × (F146 exposure by THIS observatory + overhead) = ~3h
- 1 × (F213 exposure + overhead) = ~12 min
- **Total cadence**: ~6.4 hours between F146 observations of field 0

Example Observatory Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**observatory.observatory:**

.. code-block:: text

   DIAMETER=4.0
   PIXELSCALE=0.1
   PSFFWHM=0.8
   BACKGROUND=20.0
   READOUT=3.0
   READ_OHEAD=5.0

Weather Profiles
----------------

Weather files define observing conditions over time, specifying when observations can occur.

Format
~~~~~~

**Two columns**: ``time(days) observing_flag``

- **time**: Days since simulation start
- **observing_flag**: 1 = good weather (observe), 0 = bad weather (no observation)

.. important::
   Weather files must cover the **entire simulation duration** specified by ``SIMULATION_LENGTH``
   in your parameter file. Gulls checks weather every 0.25 days (6 hours).

**Example** (``smoke.weather``):

.. code-block:: text

   0.00 1
   0.25 1
   0.50 1
   0.75 0
   1.00 1
   ...
   200.75 1

**Generating weather files**:

Use ``scripts/make_weather.py`` to generate weather profiles. Edit the script to set:

- ``init_zero``: Days of bad weather at start (burn-in period)
- ``length``: Days of good weather (observing campaign)
- ``final_zero``: Days of bad weather at end
- Sampling interval (typically 0.25 days)

Rates Files
-----------

Rates files define the expected event rate and parameter ranges for generating microlensing events.

Format
~~~~~~

**Header line** (comment): ``# field u0min u0max t0min t0max tEmin tEmax rate``

**Data lines**: One per field

+----------+--------+----------------------------------------------+
| Column   | Type   | Description                                  |
+==========+========+==============================================+
| field    | int    | Field identifier                             |
+----------+--------+----------------------------------------------+
| u0min    | float  | Minimum impact parameter                     |
+----------+--------+----------------------------------------------+
| u0max    | float  | Maximum impact parameter                     |
+----------+--------+----------------------------------------------+
| t0min    | float  | Minimum peak time (JD)                       |
+----------+--------+----------------------------------------------+
| t0max    | float  | Maximum peak time (JD)                       |
+----------+--------+----------------------------------------------+
| tEmin    | float  | Minimum Einstein crossing time (days)        |
+----------+--------+----------------------------------------------+
| tEmax    | float  | Maximum Einstein crossing time (days)        |
+----------+--------+----------------------------------------------+
| rate     | float  | Event rate (events per star per year)        |
+----------+--------+----------------------------------------------+

**Example** (``smoke.rates``):

.. code-block:: text

   # field u0min u0max t0min t0max tEmin tEmax rate
   0 0.001 1.0 2458849.0 2458949.0 10.0 100.0 1.0

.. note::
   The ``rate`` value is used as a weight for Monte Carlo sampling. Events are generated
   according to this rate and the parameter ranges specified.

Planet Files
------------

Planet files define orbital properties for planetary companions in binary-lens simulations.

Format
~~~~~~

**Header line** (comment): ``# mass(Msun) a(au) inc(deg) phase(deg)``

**Data lines**: One planet per line

+----------+--------+----------------------------------------------+
| Column   | Type   | Description                                  |
+==========+========+==============================================+
| mass     | float  | Planet mass (solar masses)                   |
+----------+--------+----------------------------------------------+
| a        | float  | Semi-major axis (AU)                         |
+----------+--------+----------------------------------------------+
| inc      | float  | Orbital inclination (degrees)                |
+----------+--------+----------------------------------------------+
| phase    | float  | Orbital phase (degrees)                      |
+----------+--------+----------------------------------------------+

**Example** (``smoke.planets.0.0``):

.. code-block:: text

   # mass(Msun) a(au) inc(deg) phase(deg)
   1.0e-6 0.5 30.0 0.0

**File naming**: ``<PLANET_ROOT>.<field>.<subrun>``

- ``PLANET_ROOT``: Specified in parameter file (e.g., ``smoke.planets``)
- ``field``: Field number
- ``subrun``: Subrun number

.. tip::
   Generate large planet catalogs efficiently using the 
   `gulls-planets <https://github.com/gulls-microlensing/gulls-planets>`_ repository,
   which supports various mass functions and orbital distributions.

Validation
----------

Gulls provides comprehensive validation tools to check your input files before running simulations. This catches common configuration errors and data format issues early, saving time and preventing failed runs.

Basic Usage
~~~~~~~~~~~

**Recommended approach** - validate using your parameter file:

.. code-block:: bash

   python3 scripts/validate_inputs.py your_parameter_file.prm

This performs comprehensive validation of all files referenced in your parameter file.

**Legacy approach** - validate individual files:

.. code-block:: bash

   python3 scripts/validate_inputs.py --sources sources.dat --lenses lenses.dat --params config.prm

**Batch validation** - check all parameter files in a directory:

.. code-block:: bash

   python3 scripts/validate_inputs.py --all

Validation Checks
~~~~~~~~~~~~~~~~~

The validation tool performs these checks in order:

1. **File Existence** - All referenced files actually exist
   - Observatory list and files
   - Source/lens/starfield catalogs and lists  
   - Weather files and directories
   - Rates files (if specified)
   - Planet files (if using planet executables)
   - Sequence files referenced in observatory files

2. **Parameter File Format** - Valid syntax and required parameters

3. **Catalog Columns** - Required column headers present
   - Standard columns: ``mul``, ``mub``, ``Mass``, ``Radius``, ``Dist``
   - Binary columns (when ``MULTIPLE_SOURCES=1``): ``Is_Binary``, ``ID``, ``primary_ID``, etc.

4. **Source/Lens Compatibility** - Valid distance relationships
   - Sources must be farther than lenses (``source.Dist > lens.Dist``)
   - Reasonable distance ranges (positive, < 50 kpc, not NaN/Inf)
   - At least some valid source/lens pairs exist

5. **NFILTERS Matching** - Parameter matches catalog format
   - ``NFILTERS`` equals number of magnitude columns in catalogs
   - Magnitude columns come before physical property columns

6. **Weather Coverage** - Weather files cover simulation duration
   - Weather file time range ≥ ``SIMULATION_LENGTH``
   - Suggests using ``scripts/make_weather.py`` if coverage insufficient

7. **Rates File Validity** - Parameter ranges are reasonable
   - ``u0min < u0max``, ``t0min < t0max``, ``tEmin < tEmax``
   - Non-negative rates and impact parameters
   - Positive Einstein times

8. **Sequence Observations** - At least one observation scheduled
   - At least one line with ``Nstack > 0`` in sequence files
   - Prevents "observatory never takes images" configuration errors

Example Output
~~~~~~~~~~~~~~

**Successful validation:**

.. code-block:: text

   Validating catalogs specified in: my_simulation.prm
   ======================================================================
   
   1. Checking parameter file format...
      ✓ Parameter file format valid
   
   2. Checking input files exist...
      ✓ All input files found
   
   3. Checking required catalog columns...
      ✓ All required columns present
   
   4. Checking source/lens compatibility...
      ✓ Valid source/lens pairs exist
      ✓ Distance values are reasonable
   
   5. Checking binary source requirements...
      ⊘ Skipped (MULTIPLE_SOURCES not enabled)
   
   6. Checking NFILTERS matches catalog format...
      ✓ NFILTERS matches magnitude columns
   
   7. Checking weather file coverage...
      ✓ Weather file covers simulation duration
   
   8. Checking rates file validity...
      ✓ Rates file parameters are valid
   
   9. Checking observing sequence has observations...
      ✓ Sequence file contains observations
   
   ======================================================================
   ✓ All validations passed!
   
   Your input files appear to be correctly formatted.

**Failed validation:**

.. code-block:: text

   ======================================================================
   ✗ Validation failed:
     Missing input files:
     - Source catalog: /path/to/nonexistent.dat
     - Observatory file: /path/to/missing.observatory
   
   Check that all file paths in your parameter file are correct.
   ======================================================================

**Warning about GSL fallbacks:**

.. code-block:: text

   ======================================================================
   ⚠ Warning: Using GSL fallback implementations for random number generation.
     The actual Numerical Recipes implementations are preferred for production runs.
     GSL fallbacks are provided for CI/testing purposes only.
   ======================================================================

Integration with CI
~~~~~~~~~~~~~~~~~~~

The same validation functions are used in:

- **Smoke tests** - Validates all test cases before running
- **CI pipeline** - Catches issues in pull requests
- **User validation** - ``scripts/validate_inputs.py`` for manual checking

This ensures consistent validation across development, testing, and production use.

Troubleshooting Validation Errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**"Missing input files"** - Check file paths in parameter file are correct

**"No valid source/lens pairs"** - Verify source distances > lens distances

**"NFILTERS mismatch"** - Count magnitude columns in your catalogs

**"Weather file coverage insufficient"** - Use ``scripts/make_weather.py`` to generate longer weather profiles

**"Sequence file has no observations"** - Add at least one line with ``Nstack > 0``

**"Binary source columns missing"** - Add required columns when ``MULTIPLE_SOURCES=1``

Examples
--------

See ``smoke_test/assets/`` for complete working examples of all input file types.
