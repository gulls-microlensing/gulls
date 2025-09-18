#include "getPlanetvals.h"
#include "buildEvent.h"
#include "dcdw.h"
#include "croin.h"
#include "structures.h"
void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
{
//extract and calculate the planet parameters

  int sdx = Event->id;
  int ln = Event->lens;
  Event->params.resize(NPLANETINPUT+NPLANETDERIV);
  Event->paramsHeader.resize(NPLANETINPUT+NPLANETDERIV);
  Event->tcroin = 0.0;
  Event->ucroin = 0.0;
  Event->rcroin = 0.0;
  //extract the input data
  for(int i=0;i<NPLANETINPUT;i++)
    {
      Event->params[i] = (*Planets)[sdx].data[i];
    }
  Event->paramsHeader[PMASS] = string("mass");
  Event->paramsHeader[AA] = string("semimajoraxis");
  Event->paramsHeader[PHASE] = string("orbphase");
  Event->paramsHeader[INC] = string("inclination");
  Event->paramsHeader[QQ] = string("q");
  Event->paramsHeader[SS] = string("s");
  Event->paramsHeader[TT] = string("period"); 

  //Semimajor axis is currently holding the a/aHZ position - swap them...
  //Event->params[AA] = Event->params[AA] * sqrt(pow(10.0,-0.4*(Lenses->data[ln][MBOL]-4.75)));

  Paramfile->parameterization=1;
  
  //Calculate the derived planet properties
  double q = Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][Lenses->MASS];

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
  double s = Event->params[SS];

  cout << Event->id << endl;
  cout << "separation: " << s << endl;
  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][Lenses->MASS]));
  //Generate a ucroin, tcroin and alphacroin and convert to u0 and alpha

  //double rcroin;
  double u0max=1;

  Event->tcroin = -1e9;
  while(!inSeason(Event->tcroin,Paramfile,World))
    Event->tcroin = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);
  Event->ucroin = u0max*(-1 + 2*ran2(Paramfile->seed));
  //cout << "check:tcroin: " << Event->tcroin << endl;
  //setupParallax(Event->tcroin, Paramfile, World, Event, Sources, Lenses);
  Paramfile->tref = Event->tcroin;

  usecroin(Event->params[SS], Event->params[QQ], Event->tcroin, Event->tE_r, Event->ucroin, Event->alpha, &Event->t0, &Event->u0, &Event->rcroin);

  //recompute u0 to deal with cases where the source is larger than rcroin
  u0max = ( Event->rs>Event->rcroin ? Event->rcroin+Event->rs : Event->rcroin)/Event->rcroin;
  Event->ucroin = u0max*(-1 + 2*ran2(Paramfile->seed));
  usecroin(Event->params[SS], Event->params[QQ], Event->tcroin, Event->tE_r, Event->ucroin, Event->alpha, &Event->t0, &Event->u0, &Event->rcroin);

  //Event->t0 = tcroin - Event->tE_r*(xcroin*cos(pi/180.0*(Event->alpha+90)) + ycroin*sin(pi/180.0*(Event->alpha+90)));
  //Event->u0 = ucroin*rcroin - (-xcroin*sin(pi/180.0*(Event->alpha+90)) + ycroin*cos(pi/180.0*(Event->alpha+90)));

  //cout << "(tcroin,ucroin)=( " << tcroin << " , " << ucroin << " )\n";
  //cout << "(t0,u0)=( " << Event->t0 << " , " << Event->u0 << " )\n";

  Event->u0max = u0max*Event->rcroin;
  Event->t0range = inSeasont0range(Paramfile,World,Event);

}
