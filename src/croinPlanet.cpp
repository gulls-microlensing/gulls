#include "getPlanetvals.h"
#include "buildEvent.h"
#include "dcdw.h"
#include "croin.h"

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

  //Semimajor axis is currently holding the a/aHZ position - swap them...
  //Event->params[AA] = Event->params[AA] * sqrt(pow(10.0,-0.4*(Lenses->data[ln][MBOL]-4.75)));

  //Calculate the derived planet properties
  double q = Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][MASS];

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
			   / (Event->params[PMASS] + Lenses->data[ln][MASS]));

  //Generate a u0croin, t0croin and alphacroin and convert to u0 and alpha

  double rcroin;
  double u0max=1;

  double t0croin = -1e9;
  while(!inSeason(t0croin,Paramfile,World))
    t0croin = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);
  double u0croin = u0max*(-1 + 2*ran2(Paramfile->seed));

  setupParallax(t0croin, Paramfile, World, Event, Sources, Lenses);

  usecroin(Event->params[SS], Event->params[QQ], t0croin, Event->tE_r, u0croin, Event->alpha, &Event->t0, &Event->u0, &rcroin);

  //Event->t0 = t0croin - Event->tE_r*(xcroin*cos(pi/180.0*(Event->alpha+90)) + ycroin*sin(pi/180.0*(Event->alpha+90)));
  //Event->u0 = u0croin*rcroin - (-xcroin*sin(pi/180.0*(Event->alpha+90)) + ycroin*cos(pi/180.0*(Event->alpha+90)));

  //cout << "(t0croin,u0croin)=( " << t0croin << " , " << u0croin << " )\n";
  //cout << "(t0,u0)=( " << Event->t0 << " , " << Event->u0 << " )\n";

  Event->u0max = u0max*rcroin;
  Event->w = Event->u0max*Event->raww / (Event->tE_r/Event->tE_h);

}
