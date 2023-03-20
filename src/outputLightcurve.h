#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "strfns.h"
#include "structures.h"

#include "info.h"
#include "mathFns.h"
#include "constdefs.h"
void outputLightcurve(struct event *Event, struct filekeywords* Paramfile, struct slcat* Sources, struct slcat* Lenses);

void outputImages(struct event *Event, struct obsfilekeywords World[], struct slcat* Sources, struct filekeywords* Paramfile);
