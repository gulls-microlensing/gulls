Running Gulls
=============

After installing Gulls and meeting all prerequisites, follow these steps to run microlensing simulations.

Quick Start
-----------

**Basic workflow:**

1. Build Gulls (see :doc:`install_gulls`)
2. Prepare input files (catalogs, observatories, parameter file)
3. Validate inputs (recommended)
4. Run simulation
5. Analyze outputs

Preparing Input Files
---------------------

Gulls requires several types of input files. See :doc:`input_formats` for detailed specifications.

Star Catalogs (Sources & Lenses)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Generate catalogs using** `Synthpop <https://synthpop.readthedocs.io/en/latest/>`_:

Synthpop is a stellar population synthesis tool for generating realistic catalogs of stars
drawn from Galactic models. For Gulls simulations:

- **Source stars**: Typically H-band magnitude limit ≤ 25 (or ≤ 27 for Roman Space Telescope fields like W146)
- **Lens stars**: No magnitude limit (mass is the key factor); typically ~10⁴ stars per field
- **Field size**: Small solid angle (δΩ ~ 10⁻⁴ deg²) representing a 0.25° × 0.25° sight line

.. important::
   Source–lens pairs where the source is closer than the lens are automatically rejected
   by Gulls. The validation tool can check this before running.

See the `Synthpop documentation <https://synthpop.readthedocs.io/en/latest/>`_ for
catalog generation instructions.

Observatory Files
~~~~~~~~~~~~~~~~~

**Pre-configured observatory files** are available in the companion repository:

`gulls-microlensing/Parameterfiles <https://github.com/gulls-microlensing/Parameterfiles>`_

This repository contains:

- Observatory configurations (telescope diameter, pixel scale, PSF, backgrounds, etc.)
- Field centres and observing sequences
- Detector properties and throughput curves
- Weather profiles

You can use these files directly or create your own following the format in :doc:`input_formats`.

Planet Files
~~~~~~~~~~~~

**Generate planet catalogs** using the companion repository:

`gulls-microlensing/gulls-planets <https://github.com/gulls-microlensing/gulls-planets>`_

Fast Python scripts that generate millions of planet entries in seconds, supporting:

- Sumi 2023 composite mass function (planets + brown dwarfs)
- Log-uniform baseline distributions
- Sensitivity-map grids (deterministic mass/semi-major-axis)
- Customizable mass ranges, orbital parameters, and inclinations

Clone the repository and follow its README for usage instructions.

Creating a Parameter File
--------------------------

The parameter file (``.prm``) orchestrates your simulation. Here's a minimal example:

.. code-block:: text

   # Basic configuration
   RUN_NAME=my_simulation
   OUTPUT_DIR=/path/to/output/
   FINAL_DIR=/path/to/final/
   EXECUTABLE=gulls_std.x

   # Input catalogs and configurations
   OBSERVATORY_DIR=observatories/
   OBSERVATORY_LIST=my_observatories.list

   SOURCE_DIR=sources/
   SOURCE_LIST=my_sources.sources
   SOURCE_COLOURS=0

   LENS_DIR=lenses/
   LENS_LIST=my_lenses.lenses
   LENS_COLOURS=0

   STARFIELD_DIR=starfields/
   STARFIELD_LIST=my_starfields.starfield

   # Optional: for planetary simulations
   PLANET_DIR=planets/
   PLANET_ROOT=planets.

   # Simulation parameters
   NFILTERS=6
   SUBRUNSIZE=1000
   NSUBRUNS=10
   PARALLAX=1
   AMIN=1.0001

See ``smoke_test/parameterfiles/`` for complete working examples.

Validating Inputs
-----------------

**Before running**, validate your catalogs to catch common issues:

.. code-block:: bash

   python3 scripts/validate_inputs.py my_simulation.prm

This checks:

- Required column headers
- Valid source/lens distance pairs (sources must be farther than lenses)
- Reasonable distance ranges (positive, < 50 kpc, not NaN/Inf)
- Binary source column requirements (when ``MULTIPLE_SOURCES=1``)
- GSL fallback detection (warns if production runs need Numerical Recipes)

Running a Simulation
--------------------

Available Executables
~~~~~~~~~~~~~~~~~~~~~

Gulls provides several executables for different simulation types:

+------------------+-------------------------------------------------------------------------+
| Executable       | Purpose                                                                 |
+==================+=========================================================================+
| gulls_std.x      | Standard planet (bound orbit) simulations with random trajectories      |
+------------------+-------------------------------------------------------------------------+
| gulls_croin.x    | Caustic Region Of INterest - planet simulations with trajectories       |
|                  | constrained to pass near caustics (speeds up planet detection studies)  |
+------------------+-------------------------------------------------------------------------+
| gullsFish.x      | Fisher matrix analysis for parameter uncertainty estimation             |
+------------------+-------------------------------------------------------------------------+

.. note::
   Additional legacy executables exist (e.g., ``gullsSingle.x``, ``gullsFFP.x``, ``gullsHZ.x``)
   but are not yet integrated into the CMake build system or smoke tests. Some may require
   updates to work with current dependencies.

Command-Line Usage
~~~~~~~~~~~~~~~~~~

Basic syntax:

.. code-block:: bash

   ./bin/gulls_std.x -i <parameter_file> -s <instance> -f <field> {-d}

**Required flags:**

- ``-i <parameter_file>`` - Path to your ``.prm`` parameter file
- ``-s <instance>`` - Simulation instance number (usually 0 for initial runs)
- ``-f <field>`` - Field number from your sources/lenses catalog

**Optional flags:**

- ``-d`` - Debug mode (repeat up to 3 times for increased verbosity: ``-d -d -d``)

**Example:**

.. code-block:: bash

   cd /path/to/gulls/
   ./bin/gulls_std.x -i parameterfiles/my_simulation.prm -s 0 -f 0 -d

This runs instance 0 of field 0 with debug output enabled.

Diagnostic Runs
~~~~~~~~~~~~~~~

For initial testing, use these settings in your parameter file:

.. code-block:: text

   PRETTY_PICS=1
   PRETTY_PICS_DIMENSIONS=128  # Or larger for better diagnostics
   OUTPUT_LC=1                  # Output every lightcurve
   OUTPUT_IMAGES=1              # Save baseline and peak images
   OUTPUT_ONALL=1              # Output all events (not just detections)

This produces maximum diagnostic output:

- **Lightcurves** (``.lc`` files): ``.det.lc`` shows detectable events, ``.all.lc`` shows all events
- **Images**: One at peak magnification, one at baseline for each filter
- Inspect images in ds9 using zscale to check for:
  
  - Realistic starfields
  - Visible microlensing events (blink between peak and baseline images for F_peak > ~1.3)

Understanding Output
--------------------

Gulls writes output to the ``OUTPUT_DIR`` specified in your parameter file:

**Key output files:**

- ``<run_name>_<field>_<instance>.out`` - Summary file with event parameters
- ``<run_name>_<field>_<instance>.lc`` - Lightcurve data
- ``<run_name>_<field>_<instance>_*.fits`` - Images (if ``OUTPUT_IMAGES=1``)

See :doc:`run_output` for detailed output specifications.

Common Parameters
-----------------

**Observing conditions:**

- ``PARALLAX`` - Enable (1) or disable (0) parallax effects
- ``LENS_LIGHT`` - Include (1) or exclude (0) lens star light (set to 0 for free-floating planets)

**Event selection:**

- ``AMIN`` - Minimum magnification threshold (e.g., 1.0001)
- ``U0MAX`` - Maximum impact parameter
- ``MIN_CHISQUARED`` - Minimum χ² for detection

**Simulation control:**

- ``SUBRUNSIZE`` - Number of events per subrun
- ``NSUBRUNS`` - Number of subruns
- ``MAXTIME`` - Maximum walltime (seconds)
- ``LC_TIMEOUT`` - Maximum time per lightcurve generation (seconds, prevents VBM hangs)

**Photometry:**

- ``IDEAL_PHOTOMETRY=-1`` - Fastest; use for initial tests
- ``LARGEPSFMAG`` - Magnitude threshold for PSF vs aperture photometry

See :doc:`parameter_reference` for comprehensive parameter documentation.

Troubleshooting
---------------

**White square images in ds9:**

1. Select "zscale" scaling in ds9
2. If still white, check for bad magnitudes in star/source/lens catalogs
3. Try larger ``PRETTY_PICS_DIMENSIONS`` (e.g., 512) to reduce impact of bright stars
4. Inspect color-magnitude diagrams of catalogs for anomalies

**"No valid fields were loaded":**

- Check that ``NFILTERS`` matches your observatory configuration
- Verify observatory files are in correct paths

**Simulation hangs or times out:**

- Check ``LC_TIMEOUT`` parameter (default 10 seconds may be too short for complex events)
- Increase for challenging binary configurations (e.g., ``LC_TIMEOUT=600``)
- Review validation warnings about catalog issues

**GSL fallback warnings:**

- For production science, replace ``src/classes/random.cpp`` and ``src/classes/zroots2.cpp``
  with licensed Numerical Recipes implementations
- GSL fallbacks are for CI/testing only

For more issues, see :doc:`basic_troubleshooting`.

Additional Resources
--------------------

- **Example configurations**: See ``smoke_test/`` for working examples
- **Postprocessing & visualization**: `gulls-postprocessing <https://github.com/gulls-microlensing/gulls-postprocessing>`_ - Jupyter notebooks for analyzing simulation output
- **Synthpop**: https://synthpop.readthedocs.io/en/latest/
- **Observatory files**: https://github.com/gulls-microlensing/Parameterfiles
- **Planet generators**: https://github.com/gulls-microlensing/gulls-planets
- **Original papers**: Penny et al. (`2013 <https://ui.adsabs.harvard.edu/abs/2013AAS...22143503P/abstract>`_, `2014 <https://ui.adsabs.harvard.edu/abs/2014ApJ...790..142P/abstract>`_, and `2019 <https://ui.adsabs.harvard.edu/abs/2019ApJS..241....3P/abstract>`_)

.. tip::
   Before postprocessing, use ``scripts/reduce_gulls.py`` to convert raw simulation output
   into HDF5 format (``.det`` and ``.out`` files) for easier analysis.
