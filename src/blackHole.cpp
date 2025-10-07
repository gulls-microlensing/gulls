#include "getPlanetvals.h"
#include "buildEvent.h"

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

  //extract the input data
  for(int i=0;i<NPLANETINPUT;i++)
    {
      Event->params[i] = (*Planets)[sdx].data[i];
    }

  //Calculate the derived planet properties
  Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][Lenses->MASS];

  double sqrtq = sqrt(Event->params[QQ]);

  Event->rE *= sqrtq;
  Event->tE_h *= sqrtq;
  Event->thE *= sqrtq;
  Event->piE /= sqrtq;
  Event->rs /= sqrtq;

  double x,y;
  x = Event->params[AA] * cos(Event->params[PHASE]*TO_RAD);
  y = Event->params[AA] * sin(Event->params[PHASE]*TO_RAD) 
    * cos(Event->params[INC]*TO_RAD);

  Event->params[SS] = sqrt(x*x + y*y) / Event->rE;

  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][Lenses->MASS]));

  //Event->t0=-1e9;
  //while(!inSeason(Event->t0,Paramfile,World))
  //  Event->t0 = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);

  //setupParallax(Event->t0, Paramfile, World, Event, Sources, Lenses);
  Paramfile->tref = Event->t0;

  //compute_u0(Paramfile, World, Event, Paramfile->seed);
  Event->u0max = ( 2*Event->rs>1 ? 2*Event->rs : 1);
  Event->u0 = (-1+2*ran2(Paramfile->seed))*Event->u0max;
  //Event->w = Event->u0max*Event->raww / (Event->tE_r/Event->tE_h); //Now dealt with by computeParallax
  Event->t0range = inSeasont0range(Paramfile,World,Event);
  Event->weight_scale = sqrt(Event->params[QQ]);
}
