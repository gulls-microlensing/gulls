#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "strfns.h"
#include "structures.h"

#include "info.h"
#include "mathFns.h"
#include "constdefs.h"
void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr);
