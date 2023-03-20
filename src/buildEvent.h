#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "structures.h"
//#include "findBlending.h"

void buildEvent(struct event *Event, struct obsfilekeywords World[], vector<vector<vector<double> > >* starfield, vector<double>* starfielddata, struct filekeywords *Paramfile, struct slcat *Sources, struct slcat *Lenses, int sdx, char* instance_, long *idum);
void addstars(struct event *Event, struct obsfilekeywords World[], 
	      vector<vector<vector<double> > >* sf, 	      
	      vector<double>* sfdata, 
	      struct filekeywords *Paramfile, struct slcat *Sources,
	      struct slcat* Lenses, long* idum);
void computeBlending(struct event *Event, struct obsfilekeywords World[],
		     struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses);
void drawsl(struct filekeywords* Paramfile,  struct obsfilekeywords World[], struct event *Event, 
	    struct slcat *Sources, struct slcat *Lenses, long* idum);
void compute_u0(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event* Event, long* idum);
int inSeason(double t0, struct filekeywords* Paramfile, struct obsfilekeywords World[]);
void setupParallax(double tref, struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event *Event, struct slcat *Sources, struct slcat *Lenses);
void setupObsGroups(struct filekeywords* Paramfile, struct event *Event);
//void buildEvent(struct event *Event, struct obsfilekeywords World[], vector<vector<vector<vector<double> > > >* starfield, vector<vector<double> >* starfielddata, struct filekeywords *Paramfile, struct galaxy *Galaxy, struct hist (*BlendData)[3], int sdx, long *idum);

