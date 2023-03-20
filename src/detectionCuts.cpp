#include "lightcurveFitter.h"
#include "detectionCuts.h"

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  Event->deterror=0;
  Event->detected=0;

  for(int obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {

      Event->currentgroup=obsgroup;

      //Determine if an event is detected 

      //Fit a PS lightcurve
      lightcurveFitter(Paramfile, Event);

      /* If finite source star effects important, fit FS lightcurve */
      if(abs(Event->umin) < 10*Event->rs && Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared)
	{
	  Event->flag_needFS[obsgroup]=1;
	  lightcurveFitter_FS(Paramfile, Event); 
	}

      //did we detect it?
      if((!Event->flag_needFS[obsgroup] && Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared)
	 || (Event->flag_needFS[obsgroup] && Event->FSPL[obsgroup].chisq>Paramfile->minChiSquared))
	{
	  Event->detected=1;
	}
    
    }
}
