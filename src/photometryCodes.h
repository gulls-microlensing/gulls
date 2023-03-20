#ifndef PHOTOMETRYCODES

//basic aperture photometry
static const int IDEAL_APERTURE = 0;
static const int APERTURE = 1;
static const int IDEAL_WEIGHTED = 2;
static const int WEIGHTED = 3;

//Fast photometry sacrifices some realism for speed - precomputing the coefficients to an error against magnitude equation, and calculating the faintest saturation which is checked, removing the need to produce images for each data point
static const int FASTAP = -1;
static const int IDEAL_FASTAP = -2;

#define PHOTOMETRYCODES
#endif
