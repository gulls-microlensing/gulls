#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "lightcurveFitter.h"
#include "integerPowers.h"

#include<fstream>

#define DEBUGVAR 0

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  double rs, Gamma;
  double amp;
  double tt, uu, u;

  Event->Amax=-1;
  Event->umin=1e50;

  int idx,obsidx;

  //work out the event parameters in the fortran parametrization

  rs = Event->rs;	/* source size */
  Gamma = Event->gamma;	/* limb-darkening profile */

  Event->lcerror=0;
  for(int obsgroup=0;obsgroup<int(Event->obsgroups.size());obsgroup++)
    {
      Event->PSPL[obsgroup].chisq=0;
      Event->flatlc[obsgroup]=0;
      for(int i=0;i<MAX_NUM_OBSERVATORIES;i++)
	Event->PSPL[obsgroup].chisqvec[i]=0;
    }

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0||Event->allsat) 
    {
      return;
    }

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];

      /* Compute the magnification */
      tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
      uu = Event->u0;

      if(Paramfile->pllxMultiplyer)
	{
	  tt += Event->pllx[obsidx].tshift(Event->jdepoch[idx]);
	  uu += Event->pllx[obsidx].ushift(Event->jdepoch[idx]);
	  //tt += Event->pllx[obsidx].tshift(idx-idxshift[obsidx]); wrong
	  //uu += Event->pllx[obsidx].ushift(idx-idxshift[obsidx]);
	}

      Event->umin=min(Event->umin,qAdd(tt,uu));

      u = qAdd(uu,tt);

      if(DEBUGVAR) cout << "witt "; cout.flush();
      amp = wittFSMagnification(u,rs);

      //cout << amp << endl;
      Event->Atrue[idx]=amp;

      //keep track of highest magnification
      if(amp>Event->Amax) 
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}

  
    }

}

