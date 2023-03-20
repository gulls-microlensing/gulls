#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include "strfns.h"
#include <time.h>

void printEpochs(struct obsfilekeywords World[], int obsidx);
void showIsseenby(struct obsfilekeywords World[], struct event *Event, int nobs);
//void writeEventHdr(struct event *Event, struct filekeywords *Paramfile,ofstream& outfile_ptr, ofstream& logfile_ptr);
void writeEventParams(struct filekeywords* Paramfile, struct event *Event, struct slcat* Sources, struct slcat* Lenses, ofstream& ofile);
void writeEventLC(struct event *Event , int idx, char eventprefix[]);
void errorHandler(int errval);
void progressbar(int idx, int niter, time_t st1);

void clock2str(time_t now, time_t start);
