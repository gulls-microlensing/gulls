Parameter Reference
===================

This document lists all parameters available in Gulls parameter files.

Basic Configuration
-------------------

RUN_NAME
~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Name for this simulation run. Used to create output directory structure.

OUTPUT_DIR
~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Base directory for output files.

FINAL_DIR
~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Final output directory (usually same as OUTPUT_DIR).

EXECUTABLE
~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Which Gulls executable to run (gulls_std.x, gulls_croin.x, gullsFish.x).

Input Directories
-----------------

OBSERVATORY_DIR
~~~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing observatory configuration files.

OBSERVATORY_LIST
~~~~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: File listing available observatories.

WEATHER_PROFILE_DIR
~~~~~~~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing weather profile files.

STARFIELD_DIR
~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing starfield files.

STARFIELD_LIST
~~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: File listing available starfields.

SOURCE_DIR
~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing source catalogs.

SOURCE_LIST
~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: File listing available source catalogs.

LENS_DIR
~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing lens catalogs.

LENS_LIST
~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: File listing available lens catalogs.

PLANET_DIR
~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Directory containing planet files.

PLANET_ROOT
~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Root name for planet files.

RATES_FILE
~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: File containing event rate information.

Output Control
--------------

OUTPUT_LC
~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Generate lightcurve files (0=no, 1=yes).

OUTPUT_IMAGES
~~~~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Generate image files (0=no, 1=yes).

PRETTY_PICS
~~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Generate pretty picture plots (0=no, 1=yes).

PRETTY_PICS_DIMENSIONS
~~~~~~~~~~~~~~~~~~~~~~
- **Type**: int
- **Default**: 32
- **Description**: Size of pretty picture plots in pixels.

OUTPUT_ONERR
~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Output files even on errors (0=no, 1=yes).

OUTPUT_ONDET
~~~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Output files only for detected events (0=no, 1=yes).

OUTPUT_ONALL
~~~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Output files for all events (0=no, 1=yes).

VERBOSITY
~~~~~~~~~
- **Type**: int
- **Default**: 4
- **Description**: Verbosity level (0=quiet, higher=more verbose).

Simulation Parameters
--------------------

SET_RANDOM_SEED_TO_CLOCK
~~~~~~~~~~~~~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Set random seed to current time (0=no, 1=yes).

RANDOM_SEED
~~~~~~~~~~~
- **Type**: int
- **Default**: 12345
- **Description**: Random number seed for reproducible results.

SIMULATION_ZERO_TIME
~~~~~~~~~~~~~~~~~~~~
- **Type**: float
- **Default**: 2458849.0
- **Description**: Reference time for simulation (JD).

SUBRUNSIZE
~~~~~~~~~~
- **Type**: int
- **Default**: 
- **Description**: Number of events to simulate.

NUM_SIM_DAYS
~~~~~~~~~~~~
- **Type**: int
- **Default**: 200
- **Description**: Number of days to simulate.

REPEAT_SEQUENCE
~~~~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Repeat observing sequence (0=no, 1=yes).

Physical Parameters
-------------------

NFILTERS
~~~~~~~~
- **Type**: int
- **Default**: 32
- **Description**: Number of filter magnitude columns in source/lens catalogs. Must match the number of magnitude columns present in your input catalogs.

AMIN
~~~~
- **Type**: float
- **Default**: 1.0001
- **Description**: Minimum impact parameter.

LARGEPSFMAG
~~~~~~~~~~~
- **Type**: float
- **Default**: 25.0
- **Description**: Magnitude limit for large PSF stars.

MIN_CHISQUARED
~~~~~~~~~~~~~~
- **Type**: float
- **Default**: 0.0
- **Description**: Minimum chi-squared for detection.

U0MAX
~~~~~
- **Type**: float
- **Default**: 1.0
- **Description**: Maximum impact parameter.

PARALLAX
~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Include parallax effects (0=no, 1=yes).

LENS_LIGHT
~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Include lens light (0=no, 1=yes).

ASTROMETRY_ON
~~~~~~~~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Enable astrometry calculations (0=no, 1=yes).

ASTROMETRIC_SYS_FLOOR
~~~~~~~~~~~~~~~~~~~~~
- **Type**: float
- **Default**: 0.0
- **Description**: Systematic astrometric error floor.

Lightcurve Generation
---------------------

LC_GEN
~~~~~~
- **Type**: int
- **Default**: 1
- **Description**: Generate lightcurves (0=no, 1=yes).

LD_GAMMA
~~~~~~~~
- **Type**: float
- **Default**: 0.00
- **Description**: Limb darkening coefficient.

VBM_RELTOL
~~~~~~~~~~
- **Type**: float
- **Default**: 1.0e-6
- **Description**: Relative tolerance for VBM calculations.

VBM_ABSTOL
~~~~~~~~~~
- **Type**: float
- **Default**: 1.0e-4
- **Description**: Absolute tolerance for VBM calculations.

LC_TIMEOUT
~~~~~~~~~~
- **Type**: float
- **Default**: 60.0
- **Description**: Timeout for lightcurve generation (seconds).

Multiplicity
------------

MULTIPLE_SOURCES
~~~~~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Enable multiple source systems (0=no, 1=yes).

MULTIPLE_LENSES
~~~~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Enable multiple lens systems (0=no, 1=yes).

Observatory Groups
------------------

OBS_GROUPS
~~~~~~~~~~
- **Type**: string
- **Default**: (ALL)
- **Description**: Observatory groups to use.

OBS_GROUP_NAMES
~~~~~~~~~~~~~~~
- **Type**: string
- **Default**: 
- **Description**: Names of observatory groups.

ERROR_SCALING
~~~~~~~~~~~~~
- **Type**: int
- **Default**: 0
- **Description**: Scale errors (0=no, 1=yes).
