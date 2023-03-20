#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "constdefs.h"
#include <gsl/gsl_matrix.h>

void jd2date(double JD, int *day, int *month, int *year, double *frac);

void lst(double JD, double EastLong, double *LST);

void eq2horiz(double incoord1, double incoord2, double JD, double Lat, double Long, char Direction, double *outcoord1, double *outcoord2);

void eq2gal(double in_coord1, double in_coord2, char dir, double *out_coord1, double *out_coord2);

void suncoo(double JD, char EquinoxType, double *RA, double *Dec, double *R,double *SL, double *EquationTime);

double obliquity(double JulianDay, char Type);

double trunc(double trunc); 

void dayofyear(int day, int month, int year, int *doy);

void refellipsoid(const char RefEllips[], double *A, double *F);

int ellipsoidchoice(const char RefEllips[]);


void geod2geoc(double longitude,  double latitude, double height, const char RefEllips[], double Geoc[], double GeocCart[]);
void mooncool(double JD , double EarthPos[] ,char Algo , double *RA, double *Dec, double *HP);

double distsp(double D1, double D2, double R1, double R2);

void moon_sky_brightness(double JD, double ObjCooRA , double ObjCooDec, double longitude, double latitude, double height, double C_Ext, double Vsky,  double *DeltaV, double *D, double *ObjMoonDist, double *K);

void mag2flux(double mag, char band, double lambda0_def, double deltalambda_def, double f_nu0_def,double *flux);

int bandcode(char band);


void gal2eclip(double l, double b, double* lambda, double* beta);
void eq2eclip(double ra, double dec, double* lambda, double* beta);
double solarlambda(double jdate);
