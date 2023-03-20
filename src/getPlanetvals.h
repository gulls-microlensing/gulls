#ifndef GETPLANETVALSHEADER

#include<vector>
#include "structures.h"
#include "constdefs.h"

void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets);

#define GETPLANETVALSHEADER
#endif
