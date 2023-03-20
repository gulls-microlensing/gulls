# gulls Documentation

# Introduction
gulls simulates large numbers of individual microlensing
events involving source and lens stars that are drawn from star
catalogs produced by a population synthesis Galactic model.
Source stars are drawn from a catalog with a faint magnitude
limit of **25 $\ge$ H vega**, and lens stars from a catalog with no
magnitude limit; source–lens pairs where the distance of the
source is less than the distance of the lens are rejected. Each
catalog is drawn from a small solid angle δΩ, but represents a
larger 0°.25 $×$ 0°.25 sight line at its specified Galactic
coordinates (ℓ, b). (Penny et al. 2019)


**For a complete background on gulls simulator please read (Penny et al. 2013) and (Penny et al. 2019).**

# Installation
1. Create environment variable with `export GULLS_BASE_DIR='/path/to/gulls/root/dir/'`
    * Can be appended to your ~/.bashrc file
2. Create an empty file with `touch ~/.gulls`
    * Not entirely sure this is necessary
3. Install [GNU Scientific Library](https://www.gnu.org/software/gnuastro/manual/html_node/GNU-Scientific-Library.html)
    * For the installation, you will have to know how many CPU threads you have
4. Install [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
    * In order for gcc to locate it, I found installing to /usr/local/ worked best; you can simply use `./configure --prefix=/usr/local`
5. Install [GFortran](https://gcc.gnu.org/wiki/GFortran)
    * If you're in Ubuntu, a slightly out-of-date version is available with `sudo apt install gfortran` and should work just fine
6. Run `./configure.sh`
7. Run `make`


# Quick Start

If you already did your research and you know how to set the parameters for gulls simulator you can just go through the following instructions and run the simulation. However, if you need more information about the parameters, We provided documentation about the parameters and sub-parameters in the next section.
  - First go to `/path/to/gulls/src/` directory and open `Makefile` using a text editor.

  ```sh
  open Makefile

  ```
  - In `Makefile`, edit the base directory `BASEDIR` and change it to the base directory path for `src` on your machine. 
  ```
  BASEDIR = /PATH/TO/gulls/src
  ```
  - Make sure the directories for `LINKERFLAGS` and `CFLAGS` are [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) and [GFortran](https://gcc.gnu.org/wiki/GFortran) path on your machine. 

  ```
  FFLAGS	= -O2 -Wall 
  CFLAGS	= -Wall -ansi -O2  -I/usr/local/Cellar/cfitsio/4.0.0/include/ -lstdc++
  #CFLAGS	= -Wall -ansi -O2
  LINKERFLAGS = -O2 -L/usr/local/Cellar/cfitsio/4.0.0/lib/ -L/usr/local/Cellar/gcc/11.2.0_3/lib/gcc/11/-lcfitsio -lm -lgfortran -lstdc++
  #LINKERFLAGS = -O2 -lcfitsio -lnsl -lm -lgfortran -lrt
  STATICFLAG = -static-libgcc -static -O2 -lstdc++
  STATICLFLAG = -O2 -lcfitsio -lgfortran -lgsl -lgslcblas -lc -lm -lpthread -lstdc++
  BINDIR = ../bin
  ```
  - If you scroll down in the `Makefile` there is a list of `Executables`:
  ```
  all: gulls gullsFFP gullshzgrid gullsFish gullsFoM gullsHZ gullsSingle gullsLimb gullscroin gullskoppa gullsffpfish gullssinglefish gullsmoa gullsBozzaOM gullsBinaryStar
  ```
  - Based on the data that you wish to generate, you might choose an `Executable` to run. Let's choose `gullsSingle` which is the executable for single lens microlensing events.
  - Now run `make gullsSingle` to create the **gullsSingle.x** executable in the `/path/to/gulls/bin` directory:
  ```sh
  make gullsSingle
  ```
  - After this step, you can clear the `make` command leftovers as it is already created the file that we need. You can do that with:
  ```sh
  make clean
  ```
  - Let's change the directory to `/path/to/gulls/bin` and run the executable:
  ```sh
  cd /path/to/gulls/bin
  ./gullsSingle.x
  ```
  - As a result, you will get a flag usage guideline:
  ```
  Usage: mabuls  -i <infile> -s <instance> {-f <field>} {-d}
  ```
  Where mabuls is the executable (gullsSingle.x) and the first flag <**infile**> is the parameter file which in our case is `PATH/singlelens.prm` (We will explain the parameter file in detail later in this documentation).<**instance**> is an integer, and we suggest to use 0 for now. <**field**> is the field number that we want to observe and for this case,we use 83. The last flag is **{-d}** which debugs while the simulator runs, therefore you do not need to change it, however, you can use it up to three times if you are generating a challenging set of data as **{-d -d -d}**.

  **If you set your parameter files correctly, gulls should start generating data in the `OUTPUT_DIR` which you set in `singlelens.prm`.**
  - Congrats! You have your very own set of data using gulls simulator!

# Parameters and Sub-parameters
Now that you know how to run gulls, let's take a couple of steps back and go through the parameters, sub-parameters, and parameter files. gulls have several inputs as parameters and sub-parameters that are included in parameter files. Including:

1. **Observatories**
2. **Weather**
3. **Starfields**
4. **Sources**
5. **Lenses**
6. **Planets**

In this section, we explain what each of these parameters represents and how you can change them to simulate your specific data set.


## Observatories

The observatories folder includes all the observatories' parameters and sub-parameter files.
* Each observatory has a `*.observatory` file which lists observatory information such as:

  * **Observatory Name**
  * **Location** -----> {LATITUDE,LONGITUDE,**SPACE** (**0 for ground-based** and **1 for space-based** observatories),ORBIT}
  * **Altitude**
  * **Telescope Data**
  * **Available Filters**
  * **Camera Data**
```
#Observatory data

#Observatory name
NAME		NRO_W

LATITUDE	99.99
LONGITUDE	99.99
SPACE		1
ORBIT		2

#Altitude above MSL 
ALTITUDE	99.99	

#Telescope data
PRIMARY		2.4
BLOCKAGE	0.3

#Available filters
FILTER		5

#Camera data
PIXELSIZE	0.11
NPIX_X		150000
NPIX_Y		150000

READ_OHEAD	5.0
NFIELDS		1

REFERENCE_TEXP	85
REFERENCE_NSTACK	1
PHOTOMETRY	-1

EXTCOEFF	0
SKY_BACKGROUND	999
```

* As well as observatory information, `*.observatory` lists a set of necessary sub-parameters such as:

  * **Field Centres** 
  * **Observing Sequence** 
  * **Detector**
  * **Throughput**
  * **Weather** 

```
#Pointing 

#Filename specifying field centres
FIELDCENTRES	nro/nro.centres

#Filename specifying observing sequence
OBSERVATION_SEQUENCE	nro/nro_W149.sequence

#Filename specifying detector parameters
DETECTOR	cycle6/cycle6_F184.detector
THROUGHPUT	cycle6/cycle6_F184.throughput

#Modifiers
WEATHER_PROFILE	WFIRST6-72.weather
```

* `*.sequence` file includes exposure times, wait times, and starfields. Here is a sample of the `sequence` file with its key:

```
#Key:
#Nstack +ve, Texp +ve  (Image being taken by this instrument)
#Nstack -ve, Texp +ve  (Image being taken by other instrument)
#Nstack -ve, Texp -ve  (No image is being taken, but something else is taking 
#	     	        time e.g. a slew or readout)


#Field	Nstack	Texp	Sum     Description
BEGIN_REPEAT	44			
0	1	46.8	2495	W149 exposure
0	-1	-862.68	2410	Other fields
END_REPEAT				
```

* `*.weather` file is a list of probabilities on given days that define if you have good or bad weather. gulls check the weather every 0.25 day (6 hours) and the **good weather** indicator is **1**, it runs the observation and if the indicator is **0**, it counts as **bad weather**, therefore, no observation would occur.  
```
0 0
0.25 0
0.5 0
0.75 1
1 1
```

## Starfields

Starfields includes catalogs of non-microlensing stars generated by **BESANCON MODEL OF STELLAR POPULATION SYNTHESIS**. Starfields are usually in four categories: 

  - **Bright** stars **H mag $\leq$ 15**
  - **Moderate 1** stars **15 $\le$ H mag $\leq$ 20**
  - **Moderate 2** stars **20 $\le$ H mag $\leq$ 25**
  - **Faint** stars **25 $\le$ H mag**


## Sources

Source stars can be drawn from starfield catalogs however they have a magnitude limit of **25 $\ge$ H vega**. In the case of W146 for the Roman Space Telescope, it can go up to **27 $\ge$ H vega**. There should be $\sim$ $10^5$ stars in the field that gulls work properly.

## Lenses

Lens stars can be drawn from starfield catalogs and they don't have any magnitude limit since mass is the main factor in microlensing events. There should be $\sim$ $10^4$ stars in the field that gulls work properly. A solid angle of $10^-4$ $deg^2$ is suggested.

## Planets

Planet parameter files carry the event parameters such as **mass, Semi major axis (a), inclination (i), and eccentricity (e)**. You can use this parameter file to generate binary lens events. 


# Parameter Files
In the quick start, we talked about `singlelens.prm` file and how essential it is to run the simulation. Here we want to learn how to set up our parameter file properly.  
  - First, start by opening the parameter file with a text editor.
  ```
  open singlelens.prm
  ```
  - Set a `RUN_NAME` and path for the `OUTPUT_DIR` and `FINAL_DIR`
  - Make sure you have the right `EXECUTABLE` in the parameter file.
  ```
  RUN_NAME= singlelens
  OUTPUT_DIR=/PATH/TO/OUTPUT_DIRECTORY/gulls_test/
  FINAL_DIR=/PATH/TO/FINAL_DIRECTORY/gulls_test/
  EXECUTABLE= gullsSingle.x
  ```
  - Now we need all the parameter files that we set up in the last part and add them to the `singlelens.prm` file:
  ```
  OBSERVATORY_DIR=observatories/
  OBSERVATORY_LIST=Farzaneh.list

  WEATHER_PROFILE_DIR=weather/

  STARFIELD_DIR=starfields/
  STARFIELD_LIST=cycle6.starlist

  SOURCE_DIR=sources/
  SOURCE_LIST=cycle6.sources
  SOURCE_COLOURS=0

  LENS_DIR=lenses/
  LENS_LIST=cycle6.lenses
  LENS_COLOURS=0

  PLANET_DIR=planets/uniform/
  PLANET_ROOT=uniform.planets.

  RATES_FILE=rates/bH.pmcorrected.rates
  
  PRINCIPLE_OBSERVATORY=0
  ```
  - When preparing a run it's a good idea to do a short run with:
  ```
  PRETTY_PICS=1
  PRETTY_PICS_DIMENSIONS=128 #Or bigger if you'd like
  OUTPUT_LC=1 #If this is less than 1, it represents a probability of output
  OUTPUT_IMAGES=1
  OUTPUT_ONALL=1
  ```
  - This should produce the maximum amount of diagnostic
    output. `OUTPUT_LC` determines the frequency with which lightcurves
    are output; `OUTPUT_ONALL` will output every event generated, but if
    set to zero, `OUTPUT_ONDET` will only output when a detection occurs
    (definition of detection depends on the type of run). `OUTPUT_IMAGES`
    saves images each time a lightcurve is an output - one at peak
    magnification, one at baseline for each filter. Inspect both the
    images and the lightcurves. .det.lc files should show some level of
    detectable microlensing event in the lightcurve, .all.lc may
    not. Images should show a starfield, and depending on the peak
    magnification, blinking between peak and base images should show the
    microlensing event (requires Fpeak>~1.3 to be easily visible).
  - If you want to generate a **PARALLAX** data, set `PARALLAX=1`, if not set it to 0.
  - You can change the `IDEAL_PHOTOMETRY` as well. `IDEAL_PHOTOMETRY=-1` gives you the fastest output. 
  - `LENS_LIGHT` can turn the lens star light off (0) and on (1). In the case of free-floating planets, you can assign `LENS_LIGHT = 0`. 
  - You can set a limit for u0, using `U0MAX`.
  - `AMIN` shows the minimum magnification.
  - There are other parameters that you can change are listed below:
  ```
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
  ```

## Checklist/Troubleshooting for regular gulls runs

- If the image is a white square, first make sure you select zscale in
**ds9**, and if that doesn't reveal stars it might indicate the presence
of bad magnitudes in the star, source, or lens catalog(s), or it's
possible that you just got unlucky and had a super-bright star land in
the image - the larger the **PRETTY_PIC**, the less likely it is that such
a star will blow out the image. Try generating a very large pretty pic
and see how many super-bright stars there are (they will look like
white squares) -- too many of these and its probably worth looking at
the color magnitude diagrams of the star lists to see that there's
nothing weird in there. 

- You can run with different levels of debug information by repeating -d flags
- If you get a message about no valid stars, check that you have
	`NFILTERS` set correctly. 
