#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "strfns.h"
#include "structures.h"
#include "constdefs.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>
#include<gsl/gsl_sf_ellint.h>
#include<gsl/gsl_mode.h>
struct parameters {double z0; double rs; double Kx; double Ky;}; 

void muVisibility(double *mu,double rs,double z0,double ld1);
int lightcurveFitter(struct filekeywords* Paramfile, struct event *Event, int enablePllx=1);

int lightcurveFitter_FS(struct filekeywords* Paramfile, struct event *Event);
double wittFSMagnification(double zsc, double rs);
