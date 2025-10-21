Input File Formats
==================

This document describes the required formats for GULLS input files.

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

Example Observatory Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**observatory.observatory:**

.. code-block:: text

   DIAMETER=4.0
   PIXELSCALE=0.1
   PSFFWHM=0.8
   BACKGROUND=20.0
   READOUT=3.0

Starfield Files
---------------

Starfield files define the distribution of background stars.

Format
~~~~~~

.. code-block:: text

   # Field ID, RA, DEC, number of stars
   0 270.0 -29.5 1000

Weather Profiles
----------------

Weather files define observing conditions over time.

Format
~~~~~~

.. code-block:: text

   # Time (days), seeing (arcsec), transparency
   0.0 1.2 0.95
   1.0 1.0 0.98
   2.0 1.5 0.90

Validation
----------

Use the validation scripts to check your input files:

.. code-block:: bash

   python3 scripts/validate_inputs.py --sources sources.dat --lenses lenses.dat

Examples
--------

See ``smoke_test/assets/`` for complete working examples of all input file types.
