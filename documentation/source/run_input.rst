Input files
===============================

.. raw:: html

   <hr>

.. _run_input_parameters:

Configuration Parameters
------------------------

The gulls simulator requires a comprehensive set of configuration parameters to accurately simulate the desired scenarios. The configuration is facilitated through a YAML file format, offering flexibility and ease of use for specifying simulation parameters.


.. admonition:: Configuration Parameters Overview

   The configuration file is divided into several sections, each corresponding to different components of the simulation environment:

   * :ref:`Observation <run_input_parameters_observatories>`   
   * Star Catalogs
   * Planets

For an in-depth explanation of these parameters and how they influence the simulation, consult the documentation provided in the following sections. Further details on specific simulation mechanics, including reference frames, PSF radial dependency, and platform jitter, can be found in the gulls technical documentation and relevant scientific papers.

.. raw:: html

   <hr>
.. _run_input_parameters_observatories:

*Observatories*
...............

The gulls simulator employs a detailed configuration for observatories, tailored to simulate observations accurately. Below is an example configuration for a sample observatory.

.. code-block:: shell

   #INSTRUMENT - WFIRST
   #FILTER - W169

   #Observatory data

   #Observatory name
   NAME        NRO_W

   LATITUDE    99.99
   LONGITUDE   99.99
   SPACE       1  # Allowed values: 0 (ground-based), 1 (space-based)
   ORBIT       5  # Allowed values: integers representing orbit types

   #Altitude above MSL 
   ALTITUDE    99.99    

   #Telescope data
   PRIMARY     2.4
   BLOCKAGE    0.3

   #Available filters
   FILTER      0  # Allowed values: integers corresponding to specific filters

   #Camera data
   PIXELSIZE   0.11
   NPIX_X      550000
   NPIX_Y      550000

   READ_OHEAD  5.0
   NFIELDS     1

   REFERENCE_TEXP      85
   REFERENCE_NSTACK    1
   PHOTOMETRY          -1

   EXTCOEFF    0
   SKY_BACKGROUND 999

   #Pointing 

   #Filename specifying field centres
   FIELDCENTRES nro/nro.centres

   #Filename specifying observing sequence
   OBSERVATION_SEQUENCE cycle7/c7_Z087_47s.sequence

   #Filename specifying detector parameters
   DETECTOR    cycle7/cycle6_Z087.detector
   THROUGHPUT  cycle7/cycle6_Z087.throughput

   #Modifiers
   WEATHER_PROFILE WFIRST6-72.weather

**Note on Configuration Parameters:**

This configuration exemplifies setting up an observatory within the gulls simulator, demonstrating how various aspects of the observatory and observational instruments are defined. Users are encouraged to modify these parameters based on their specific simulation requirements, ensuring the setup accurately reflects the intended observational scenario.

- **SPACE**: :guilabel:`Allowed values`

 Indicates whether the observatory is ground-based or space-based. It accepts `0` for ground-based observatories and `1` for space-based platforms.

- **ORBIT**: :guilabel:`Allowed values`

Represents the specific orbit type for space-based observatories. It accepts integer values, with each integer signifying a different orbit classification.
.. code-block:: text

   case 0: earth
   case 1: Geosynchronous
   case 2: L2
   case 3: Mars
   case 4: Jupiter
   case 5: JWST

- **FILTER**: Corresponds to the specific filter used during the observation. Integer values are used to denote to a specific filter's column number in the star catalogs.

**Location and Pointing**:

Observatory location, defined by LATITUDE and LONGITUDE, alongside the SPACE and ORBIT parameters, sets the observational platform's position. POINTING parameters, including FIELDCENTRES and OBSERVATION_SEQUENCE, DETECTOR, THROUGHPUT, and WEATHER_PROFILE are other configration sub-parameters.

### Dependency Files
.. raw:: html

   <hr>

.. _run_input_parameters_detector:

Detector
.........

The detector configuration plays a critical role in accurately simulating observational data. Below is a detailed breakdown of the detector settings for a sample instrument. These settings are specified in the associated detector file and reflect key characteristics such as bias, readout noise, thermal flux, and pixel scale.

.. code-block:: shell

   # Instrument and Band information
   Instrument: WFIRST
   Band: W149

   # Detector settings based on WFIRST reference information
   BIAS:            1000          # Bias level in counts per pixel
   READOUT:         12.12        # Read out noise in counts per pixel (Spec)
   THERMAL:         0.030        # Thermal flux in counts per pixel per sec (284K)
   DARKCURRENT:     0.1          # Dark current in counts per pixel per sec (Spec)

   PIXELSCALE:      0.11         # Pixel size in arcsec
   PSFFWHM:         0.2          # PSF FWHM in arcsec
   PSFFILE:         psfs/Cycle5_SCA09_Z087_M2V_psf.psf
   PSFSCALE:        0.0122222222 # Spacing between samples of the numerical PSF
   KERNSIZE:        84           # PSF Kernel size
   SUBPIX:          9            # Number of sub-pixel points used to place stars

   SYSTEMATIC:      0.001        # Systematic photometry error
   BITDEPTH:        679000       # Number of bits per pixel
   FULLWELL:        679000       # Number of electrons a pixel can hold
   GAIN:            1            # Inverse gain
   BLEEDING:        0            # No charge bleeding
   BLEEDACROSS:     0            # Nothing crosses columns

   DIAMETER:        2.36         # Telescope diameter in metres
   BLOCKAGE:        0.3          # Fractional linear blockage

   ZEROMAG:         26.3866      # Magnitude at which zeropoint is defined
   ZEROFLUX:        1            # Photons per second from a reference source

   BACKGROUND:      999          # Background magnitude in mags per sq arcsec
   APERTURE:        0.2          # Aperture radius (not diameter) in arcsec CHECK

   PIXELSIZE:       18           # Pixel size in microns
   CRFLUX:          0            # Cosmic ray flux in hits m^-2 s^-1 CHECK

### Note:

- The configuration parameters provided in the detector file are instrumental in defining the optical and electronic characteristics of the simulated observation. They directly influence the quality and accuracy of the simulated data.
- It is imperative to ensure that the values provided in the `PSFFILE`, `BACKGROUND`, and other parameters are aligned with the specific requirements of the simulation scenario being modeled.

.. raw:: html

   <hr>

.. _run_input_parameters_center:

Field centers
..............

Depending on how many fields telescope has is planning to observe the center of each field is defined in each lines, here is an example for "nro.centres":


.. code-block:: shell

   1.1	-1.7

.. raw:: html

   <hr>

.. _run_input_parameters_sequence:

Observing sequence
....................

The gulls simulator relies on specific sequence files to detail the observing schedule, including exposure times, stack numbers, and descriptions. The followwing `.sequence` file, for instance, outlines the observing sequence for the 7 field in 2 filters.

.. code-block:: shell

   # Observing sequence for the 7 field survey
   
   
   # Key:
   # Nstack +ve, Texp +ve  (Image being taken by this instrument)
   # Nstack -ve, Texp +ve  (Image being taken by other instrument)
   # Nstack -ve, Texp -ve  (No image is being taken, but something else is taking 
   #                        time e.g. a slew or readout)
   
   # Additional time is inserted between images of stacks by exigere for readout and dithering. 
   # The readout time for this can be set in the .observatory observatory files.
   # All other readout, slewing, filter changes etc must be specified here.
   
   # Field  Nstack  Texp    Sum     Description
   
   # F087 Set
   0       -1      286.0   366     F087 exposure
   0       -1      -2297.88 76     Other fields
   
   # W149/W169 set
   BEGIN_REPEAT    44
   0       1       46.8    2495    W149 exposure
   0       -1      -862.68  2410    Other fields
   END_REPEAT

### Notes:

- **Sequence File Customization**: The `.sequence` file demonstrates how to configure detailed observing sequences. Adjustments may be required to align with specific simulation goals, available instruments, and observational constraints.

.. raw:: html

   <hr>

.. _run_input_parameters_weather:

Weather
........

.. code-block:: shell

   Weather:
     Profile: WFIRST6-72.weather

Weather conditions play a crucial role in observational simulations. The **Weather** section allows you to specify if the telescope observes based on the weather. Each Line is a quart of day.
.. code-block:: shell

   0 0
   0.25 0
   0.5 0
   0.75 1
   1 1

**Weather**: :guilabel:`Allowed values`
It accepts `0` for NOT OBSERVING and `1` for OBSERVING.

.. raw:: html

   <hr>

.. _run_input_stars:

*Star catalogs*
----------------

gulls requires three catagory of star catalogs to draw from:

* Source catalog
* Lense catalog
* Background starfield catalogs

.. raw:: html

   <hr>


*Source catalog*
................

Source stars can be drawn from source catalogs generated by synthetic population models. Depending on users simulation goals, source stars can have a magnitude limit. For instance, in the case of W146 for the Roman Space Telescope, source's magnitude limit can go up to **27 $\ge$ H vega**. There should be $\sim$ $10^5$ stars in the observing field that gulls work properly.

.. raw:: html

   <hr>


*Lens catalog*
..............
Lens stars can be drawn from lens catalogs  generated by synthetic population models. Lens stars don't have any magnitude limit since mass is the main factor in microlensing events. There should be $\sim$ $10^4$ stars in the field that gulls work properly. A solid angle of $10^-4$ $deg^2$ is suggested.

.. raw:: html

   <hr>


*Starfield catalog*
...................

Starfield catalogs include non-microlensing stars generated by **BESANCON MODEL OF STELLAR POPULATION SYNTHESIS**. Starfields are usually in four categories: 

  - **Bright** stars **H mag $\leq$ 15**
  - **Moderate 1** stars **15 $\le$ H mag $\leq$ 20**
  - **Moderate 2** stars **20 $\le$ H mag $\leq$ 25**
  - **Faint** stars **25 $\le$ H mag**


.. raw:: html

   <hr>

.. _run_input_planets:

*Planets*
---------

Planet files carry the event parameters such as **mass, Semi major axis (a), inclination (i), and eccentricity (e)**. You can use planet file to generate binary lens events.

 .. warning::
   Planet files dictate the quantity of events that gulls simulates, including scenarios where no planet injection occurs within the simulation... raw:: html

.. raw:: html

   <hr>

.. _run_input_files:

*Parameter File*
------------------

The gulls Simulator relies on parameter files to configure simulations accurately. The `.prm` file serves as a crucial input, containing essential settings that define how your simulation operates. This document provides guidance on configuring your `.prm` effectively.

.. contents::
   :local:
   :depth: 2


1. **Opening the Parameter File**:

   Begin by locating and opening your `.prm` file in a text editor. This file dictates key simulation parameters and paths.

   .. code-block:: shell

      open singlelens.prm

2. **Setting Basic Information**:

   The parameter file requires you to specify several foundational settings:

   - **RUN_NAME**: A unique identifier for your simulation run.
   - **OUTPUT_DIR**: Directory where simulation outputs will be stored.
   - **FINAL_DIR**: Directory for storing final simulation results.
   - **EXECUTABLE**: The specific gulls executable to be used for this simulation.

   Example configuration:

   .. code-block:: shell

      RUN_NAME= singlelens
      OUTPUT_DIR=/PATH/TO/OUTPUT_DIRECTORY/gulls_test/
      FINAL_DIR=/PATH/TO/FINAL_DIRECTORY/gulls_test/
      EXECUTABLE= gullsSingle.x

3. **Specifying Simulation Components**:

   Include paths to directories and lists for various simulation components:

   - Observatories, weather conditions, starfields, sources, lenses, and planet parameters are specified through directory paths and list names.

   Example paths:

   .. code-block:: text

      OBSERVATORY_DIR=observatories/
      OBSERVATORY_LIST=observatory.list
      WEATHER_PROFILE_DIR=weather/
      STARFIELD_DIR=starfields/
      STARFIELD_LIST=cycle6.starlist
      SOURCE_DIR=sources/
      SOURCE_LIST=cycle6.sources
      LENS_DIR=lenses/
      LENS_LIST=cycle6.lenses
      PLANET_DIR=planets/uniform/
      PLANET_ROOT=uniform.planets.
      RATES_FILE=rates/bH.pmcorrected.rates
      PRINCIPLE_OBSERVATORY=0

4. **Optimizing Output for Analysis**:

   To enhance the analysis of simulation results, consider enabling additional output options:

   - **PRETTY_PICS**: Generates detailed graphical representations.
   - **OUTPUT_LC**: Controls the frequency of lightcurve outputs.
   - **OUTPUT_ONALL**: Controls the frequency of all outputs, even if there is no detection.
   - **OUTPUT_DET**: Controls the frequency of microlensing detection outputs
   - **OUTPUT_IMAGES**: Saves images at key simulation points.

   Diagnostic settings:

   .. code-block:: text

      PRETTY_PICS=1
      PRETTY_PICS_DIMENSIONS=128
      OUTPUT_LC=1
      OUTPUT_IMAGES=1
      OUTPUT_DET= 1
      OUTPUT_ONALL=1


5. **Advanced Configuration Options**:

   For specialized simulations, you might adjust the following:

   - **IDEAL_PHOTOMETRY**, **PARALLAX**, and **LENS_LIGHT** for specific observational effects.
   - **RANDOM_SEED** and timing settings to influence the stochastic elements of the simulation.

   Example advanced settings:

   .. code-block:: text

      IDEAL_PHOTOMETRY=-1
      PARALLAX=1
      LENS_LIGHT=1
      SET_RANDOM_SEED_TO_CLOCK=1
      RANDOM_SEED=1001
      SIMULATION_ZERO_TIME=2458234.000000
      SIMULATION_LENGTH=5.50803
      SUBRUNSIZE=1000
      NSUBRUNS=40
      MAXTIME=14000.000000
      NFILTERS=6
      AMIN=1.0000001
      LARGEPSFMAG=19.0
      MIN_CHISQUARED=60
      U0MAX=1.0

.. note::

   The correct configuration of the `.prm` file is vital for the successful execution of gulls simulations. Ensure all paths and settings accurately reflect your simulation goals and system setup.
