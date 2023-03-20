#include "lightcurveFitter.h"

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  //Determine if an event is detected 

  //Fit a PS lightcurve
  Event->deterror=0;
  Event->detected=0;

  if(Event->flatchi2>500)
    {
      Event->detected=1;
    }

  //Event->flatchi2 = Event->flatlc;
}
