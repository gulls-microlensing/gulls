#include "getPlanetvals.h"
#include "buildEvent.h"
#include "croin.h"

void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
{
  //extract and calculate the planet parameters

  int sdx = Event->id;
  int ln = Event->lens;
  Event->params.resize(NPLANETINPUT+NPLANETDERIV);
  Event->paramsHeader.resize(NPLANETINPUT+NPLANETDERIV);

  Event->paramsHeader[PMASS] = string("mass");
  Event->paramsHeader[AA] = string("semimajoraxis");
  Event->paramsHeader[PHASE] = string("orbphase");
  Event->paramsHeader[INC] = string("inclination");
  Event->paramsHeader[QQ] = string("q");
  Event->paramsHeader[SS] = string("s");
  Event->paramsHeader[TT] = string("period"); 

  Paramfile->parameterization=0;
  
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

  //setupParallax(Event->t0, Paramfile, World, Event, Sources, Lenses);
  Paramfile->tref=Event->t0;
  Paramfile->parameterization=0;

  //Event->w /= (Event->tE_r/Event->tE_h);

  //The croin parameters will not be accurate if there is orbital motion or parallax
  croinparam(Event->params[SS], Event->params[QQ], Event->t0, Event->tE_r, Event->u0, Event->alpha, &Event->tcroin, &Event->ucroin, &Event->rcroin);

//Place some conditions to speed things up...

  //double alp = (Event->alpha>=180?360-Event->alpha:Event->alpha);

  /*if(Event->u0>0.02)
    {
      if(Event->params[SS]<0.1) Event->allsat=1;
      if(Event->params[SS]>7.0 && Event->alpha>100 && Event->alpha<250 && abs(Event->u0)>0.5) Event->allsat=1;
      if(log10(Event->params[QQ])<-4.5)
	{
	  if(alp>90.5&&(Event->params[SS]>1.3||Event->params[SS]<0.24)) Event->allsat=1;
	  if(alp<89.5&&Event->params[SS]<0.75) Event->allsat=1;
	  if(Event->alpha<89.5 && Event->params[SS] > 3.05 + 1.0/((90-alp)/55.0)) Event->allsat=1;
	}
	}*/

}
