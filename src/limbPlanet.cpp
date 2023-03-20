#include "getPlanetvals.h"
#include "buildEvent.h"

void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
{
  //extract and calculate the planet parameters

  int sdx = Event->id;
  int ln = Event->lens;
  Event->params.resize(NPLANETINPUT+NPLANETDERIV);

  //extract the input data
  for(int i=0;i<NPLANETINPUT;i++)
    {
      Event->params[i] = (*Planets)[sdx].data[i];
    }

  //Calculate the derived planet properties
  Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][MASS];

  double sqrt1pq = sqrt(1+Event->params[QQ]);

  Event->rE *= sqrt1pq;
  Event->tE_h *= sqrt1pq;
  Event->thE *= sqrt1pq;
  Event->piE /= sqrt1pq;
  Event->rs /= sqrt1pq;
  //do not update the event rate weighting - we only update the quantities 
  //where the calculation of the lightcurve requires a normalized mass

  double x,y;
  x = Event->params[AA] * cos(Event->params[PHASE]*TO_RAD);
  y = Event->params[AA] * sin(Event->params[PHASE]*TO_RAD) 
    * cos(Event->params[INC]*TO_RAD);

  Event->params[SS] = sqrt(x*x + y*y) / Event->rE;

  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][MASS]));


  if(Event->u0max<0.01)
    {
      Event->u0=1;
      Event->w=0;
    }
  else
    {
      while(abs(Event->u0)<0.01)
	{
	  Event->u0 = Event->u0max*(2*ran2(Paramfile->seed)-1);
	}
      Event->w = (Event->u0max-0.01)*Event->raww;
    }

  setupParallax(Event->t0, Paramfile, World, Event, Sources, Lenses);
}
