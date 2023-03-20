#include "getPlanetvals.h"
#include "buildEvent.h"

void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
{
  //extract and calculate the planet parameters

  //single lenses - ignore planets, but still calculate the planet parameters to fill out the output

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

  double x,y;
  x = Event->params[AA] * cos(Event->params[PHASE]*TO_RAD);
  y = Event->params[AA] * sin(Event->params[PHASE]*TO_RAD) 
    * cos(Event->params[INC]*TO_RAD);

  Event->params[SS] = sqrt(x*x + y*y) / Event->rE;

  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][MASS]));

  Event->t0=-1e9;
  while(!inSeason(Event->t0,Paramfile,World))
    Event->t0 = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);

  //Overwrite the work of compute_u0

  double umaxmax=1.0;
  double umaxmin=0.01;

  int obsidx = Paramfile->principle_observatory;
  int filter = World[obsidx].filter;

  //Calculating u0

  //interpret Amin as the faintest primary-band peak magnitude for selection
  double Imax = Paramfile->Amin;
  double Ibase;

  //Compute using the total blending in primary band (lens+blend)

  double fs = Event->fs[obsidx];
  Ibase = Sources->mags[Event->source][filter] + 2.5*log10(fs);
  double mumin, u0max;

  if(Ibase<Imax)
    {
      u0max=1;
    }
  else
    {
      mumin = (pow(10, -0.4*(Imax-Ibase)) - 1+fs)/fs;      //(Amin-1+fs)/fs;

      //ensure that the microlensing event can be seen
      u0max = sqrt(2.0*sqrt(1.0+1.0/(mumin*mumin-1.0)) - 2.0);
      //u0max = (u0max>umaxmax?umaxmax:u0max);
      //u0max = (u0max<umaxmin?umaxmin:u0max);
    }

  Event->u0 = u0max*ran2(Paramfile->seed);

  Event->w = u0max*Event->raww;
  Event->u0max = u0max;

  setupParallax(Event->t0, Paramfile, World, Event, Sources, Lenses);
  Event->w /= (Event->tE_r/Event->tE_h);

  //Compute some interesting quantities and store them in the planet variables
  
  //time with I<15.5
  Event->params[TT] = 2 * u0max * Event->tE_r;

}


