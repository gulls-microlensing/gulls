# gulls

A microlensing simulator optimized for space-based microlensing surveys.

Warning: This repository is not yet fully operational due to our
efforts to obey licencing requirements. Most of the gulls code is
here, and if you can't wait ~week for us to complete a review to
replace code with licensing that is too restrictive, please contact
penny1@lsu.edu
   
 

## Installation Instructions

1. Create environment variable with `export GULLS_BASE_DIR='/path/to/gulls/root/dir/'`
    * Can be appended to your ~/.bashrc file
2. Create empty file with `touch ~/.gulls`
    * Not entirely sure this is necessary
3. Install [GNU Scientific Library](https://www.gnu.org/software/gnuastro/manual/html_node/GNU-Scientific-Library.html)
    * For the installation, you will have to know how many CPU threads you have
4. Install [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
    * In order for gcc to locate it, I found installing to /usr/local/ worked best; you can simply use `./configure --prefix=/usr/local`
5. Install [GFortran](https://gcc.gnu.org/wiki/GFortran)
    * If you're in Ubuntu, a slightly out of date version is available with `sudo apt install gfortran` and should work just fine
6. Run `./configure.sh`
7. Run `make`


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

## Course Files for Image Simulations

To improve separation ability, dithering may be used by `freeColour` to remove noise and artifacts. This is enabled by custom `.course` files which list sets of day,x,y,theta waypoints for images to be captured. These `.course` files are specified in each line of the detector list file.

Day is simply the timestamp of the image capture in days; the baseline day 0 is arbitrary. Lateral x,y coordinates are doubles in terms of pixel coordinates. For example, 0.5 represents a half pixel offset from the origin. Angle theta is measured as the CCW (counter-clockwise) angle of the images axes relative to the x,y grid in degrees. Rotation is centered around the exact center of the image.

The day,x,y,theta coordinates are separated by a space and each line ends in a new line character. A list of these is used to specify the sequential path the dithering follows. For example, one such path could be

```
32.5 0.25 0.25 30.0 
33.6 0.6 0.0 19.2 
35.4 -1.1 0.9 51.4 
```

Note how the (x,y) origin doesn't have to be the first point, or even used at all. Similarly, the theta doesn't even have to be set to 0 degrees, nor does the day have to start at 0 days. Since offsets are relative and star layouts are random, absolute coordinates are meaningless. Also, even negative numbers can be used.

NOTE: A `.course` file MUST be supplied, even if it is trivial. Here is an example of a trivial path for stationary images

```
0.0 0.0 0.0 0.0 
```
