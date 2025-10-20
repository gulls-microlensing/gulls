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

Required Columns
~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| RA2000.0       | float  | Right ascension (J2000)        | degrees  |
+----------------+--------+--------------------------------+----------+
| DEC2000.0      | float  | Declination (J2000)            | degrees  |
+----------------+--------+--------------------------------+----------+
| Dist           | float  | Distance to source             | kpc      |
+----------------+--------+--------------------------------+----------+
| Mass           | float  | Stellar mass                   | M☉       |
+----------------+--------+--------------------------------+----------+
| Radius         | float  | Stellar radius                 | R☉       |
+----------------+--------+--------------------------------+----------+

Optional Columns
~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| mul            | float  | Proper motion in l             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| mub            | float  | Proper motion in b             | mas/yr   |
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

Example Source Catalog
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   RA2000.0 DEC2000.0 Dist Mass Radius mul mub Vr [Fe/H] Teff logg Is_Binary ID primary_ID combined_logP
   270.1234 -29.5678 8.5 1.0 1.0 0.0 0.0 0.0 0.0 5778 4.44 0 1 0 0.0
   270.1235 -29.5679 8.5 0.8 0.8 0.0 0.0 0.0 0.0 5000 4.5 1 2 0 1.5
   270.1236 -29.5680 8.5 0.6 0.6 0.0 0.0 0.0 0.0 4500 4.6 2 3 2 1.5

Lens Catalogs
-------------

Lens catalogs contain information about potential lensing objects.

Required Columns
~~~~~~~~~~~~~~~~

+----------------+--------+--------------------------------+----------+
| Column         | Type   | Description                    | Units    |
+================+========+================================+==========+
| RA2000.0       | float  | Right ascension (J2000)        | degrees  |
+----------------+--------+--------------------------------+----------+
| DEC2000.0      | float  | Declination (J2000)            | degrees  |
+----------------+--------+--------------------------------+----------+
| Dist           | float  | Distance to lens               | kpc      |
+----------------+--------+--------------------------------+----------+
| Mass           | float  | Lens mass                      | M☉       |
+----------------+--------+--------------------------------+----------+
| mul            | float  | Proper motion in l             | mas/yr   |
+----------------+--------+--------------------------------+----------+
| mub            | float  | Proper motion in b             | mas/yr   |
+----------------+--------+--------------------------------+----------+

Example Lens Catalog
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   RA2000.0 DEC2000.0 Dist Mass mul mub
   270.1200 -29.5600 4.0 0.3 2.0 1.0
   270.1300 -29.5700 6.0 0.5 1.5 0.8

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
