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
  double q = pow(10,Event->params[PMASS]);
  Event->params[QQ] = q;
  double s = pow(10,Event->params[AA]);
  Event->params[SS] = s;

  double sqrt1pq = sqrt(1+Event->params[QQ]);
  
  cout << Event->rE << " " << q << " " <<  Event->params[QQ] << " " << sqrt1pq << endl;
    
  Event->rE *= sqrt1pq;
  Event->tE_h *= sqrt1pq;
  Event->thE *= sqrt1pq;
  Event->piE /= sqrt1pq;
  Event->rs /= sqrt1pq;
  //do not update the event rate weighting - we only update the quantities 
  //where the calculation of the lightcurve requires a normalized mass
  
  Event->params[PMASS] = q * Lenses->data[ln][MASS];
  Event->params[AA] = Event->params[SS] * Event->rE;

  double x,y;
  x = Event->params[SS] / cos(-Event->params[PHASE]*TO_RAD);
  y = Event->params[SS] / sin(-Event->params[PHASE]*TO_RAD) 
    / cos(Event->params[INC]*TO_RAD);

  Event->params[AA] = sqrt(x*x + y*y) * Event->rE;

  cout << Event->id << endl;

  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][MASS]));

  //Draw a random u0,t0, then see if it is in the croin regime
  //Generate a u0croin, t0croin and alphacroin and convert to u0 and alpha

  double rcroin,uc,tc;
  double u0max=Paramfile->u0max;

  
  double t0 = -1e9;
  while(t0<=0 || (t0>549.75+112.5 && t0<1388.0-112.5) || t0>1825.25+112.5)
    t0 = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);
  double u0 = u0max*(-1 + 2*ran2(Paramfile->seed));

  Event->t0 = t0;
  Event->u0 = u0;

  setupParallax(t0, Paramfile, World, Event, Sources, Lenses);

  croinparam(Event->params[SS], Event->params[QQ], t0, Event->tE_r, u0, Event->alpha, &tc, &uc, &rcroin);

  Event->params.push_back(uc);
  Event->params.push_back(tc);
  Event->params.push_back(inSeason(tc,Paramfile,World));

  cout << Event->params[SS] << " " << Event->params[QQ] << " " << Event->params[PMASS] << " " << Event->params[AA] << " " << t0 << " " << u0 << " " << Event->params[PHASE] << " " << Event->params[INC] << " " << Event->rE << endl;
  cout << Event->id << " " << Event->instance << " " << Event->field << " " << uc << " " << tc << " " << " " << uc*rcroin << " " << inSeason(tc,Paramfile,World) << endl;

  //Event->t0 = t0croin - Event->tE_r*(xcroin*cos(pi/180.0*(Event->alpha+90)) + ycroin*sin(pi/180.0*(Event->alpha+90)));
  //Event->u0 = u0croin*rcroin - (-xcroin*sin(pi/180.0*(Event->alpha+90)) + ycroin*cos(pi/180.0*(Event->alpha+90)));

  //cout << "(t0croin,u0croin)=( " << t0croin << " , " << u0croin << " )\n";
  //cout << "(t0,u0)=( " << Event->t0 << " , " << Event->u0 << " )\n";

  //Event->u0max = u0max*rcroin;
  Event->w = Event->u0max*Event->raww / (Event->tE_r/Event->tE_h);

}
