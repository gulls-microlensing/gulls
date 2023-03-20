#include "getPlanetvals.h"
#include "buildEvent.h"
#include "hz.h"
#include "dcdw.h"
#include "croin.h"

double poly(double x, int n, double* coeff)
{
  double result=coeff[n];
  for(int i=n-1;i>=0;i--)
    {
      result=x*result+coeff[i];
    }
  return result;
}

void getPlanetvals(struct event* Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses, vector<struct pcat>* Planets)
{
//extract and calculate the planet parameters

  int sdx = Event->id;
  int ln = Event->lens;
  int sn=Event->source;
  int obsidx;
  Event->params.resize(NPLANETINPUT+NPLANETDERIV);

  double h[11]={8.88832533333353,29.2789080330271,-444.498144552323,
		2386.58819281151,-7025.89991134296,12444.5766946781,
		-13808.4580644194,9666.9810224298,-4146.99868105641,
		995.47659849976,-102.430192057616}; //H polynomial coefficients
  double I[11]={13.1862177834363,-15.8752142623272,-65.8470509443578,
		622.500039103007,-2097.44301285205,3820.57789955284,
		-4187.30617436531,2848.24698898564,-1179.61240696304,
		272.861018214196,-27.0608201379289}; //I polynomial coefficients
  double lt[11]={3.56975498151499,-1.63133182995477,15.3300830495442,
		 -64.8030182373953,152.045848973286,-210.317900750694,
		 178.430120721367,-93.6045629118416,29.4347318019326,
		 -5.03580513420022,0.354508460472279}; //Teff polynomial coefficients
  double lb[11]={-3.53498337086972,6.33036302190073,25.6270696107729,
		 -254.499055767534,886.630620668837,-1671.18548690424,
		 1899.06103024211,-1340.78455994293,576.428087819269,
		 -138.330694295787,14.2181271245544}; //lbol polynomial coefficients

  //extract the input data
  for(int i=0;i<NPLANETINPUT;i++)
    {
      Event->params[i] = (*Planets)[sdx].data[i];
    }

  //Convert mass
  double mindex = Event->params[PMASS];
  Event->params[PMASS] = 3.00374072e-6*pow(10,0.5*int(mindex/50.0));
  Lenses->data[ln][MASS] = mindex - 50*int(mindex/50.0);
  Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][MASS];

  //Replace the lens star on the image
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    { 
      int filter = World[obsidx].filter;
      if(Paramfile->lenslight)
	{
	  World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   Lenses->mags[Event->lens][filter]);
	}
    }

  Lenses->mags[ln][3] = poly(Lenses->data[ln][MASS],10,h) + 5*log10(Lenses->data[ln][DIST]/0.01) + Lenses->data[ln][AV]*0.184 + 1.39;
  Lenses->mags[ln][1] = poly(Lenses->data[ln][MASS],10,I) + 5*log10(Lenses->data[ln][DIST]/0.01) + Lenses->data[ln][AV]*0.527 + 0.45;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    { 
      int filter = World[obsidx].filter;
      if(Paramfile->lenslight)
	{
	  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   Lenses->mags[Event->lens][filter]);
	}
    }
  computeBlending(Event, World, Paramfile, Sources, Lenses);

  double u0max=1;

  double sqrt1pq = sqrt(1+Event->params[QQ]);

  double rEold=Event->rE;
  Event->rE = rEsun*sqrt(Lenses->data[ln][MASS]*Lenses->data[ln][DIST]*(1-Lenses->data[ln][DIST]/Sources->data[sn][DIST]))*sqrt1pq;
  Event->tE_h = Event->tE_h*(Event->rE/rEold);
  Event->thE = Event->thE*(Event->rE/rEold);
  Event->piE = Event->piE/(Event->rE/rEold); 
  Event->rs = Event->rs/(Event->rE/rEold);
  Event->raww = Event->raww*(Event->rE/rEold);
  //do not update the event rate weighting - we only update the quantities 
  //where the calculation of the lightcurve requires a normalized mass

  double x,y;
  x = Event->params[AA] * cos(Event->params[PHASE]*TO_RAD);
  y = Event->params[AA] * sin(Event->params[PHASE]*TO_RAD) 
    * cos(Event->params[INC]*TO_RAD);

  Event->params[SS] = sqrt(x*x + y*y) / Event->rE;

  Event->params[TT] = sqrt(cube(Event->params[AA]) 
			   / (Event->params[PMASS] + Lenses->data[ln][MASS]));


  //Generate a u0croin, t0croin and alphacroin and convert to u0 and alpha

  double t0croin = -1e9;
  while(!inSeason(t0croin,Paramfile,World))
    t0croin = double(Paramfile->NUM_SIM_DAYS)*ran2(Paramfile->seed);
  double u0croin = u0max*(-1 + 2*ran2(Paramfile->seed));
  double rcroin;

  setupParallax(t0croin, Paramfile, World, Event, Sources, Lenses);

  usecroin(Event->params[SS], Event->params[QQ], t0croin, Event->tE_r, u0croin, Event->alpha, &Event->t0, &Event->u0, &rcroin);

  //Event->t0 = t0croin - Event->tE_r*(xcroin*cos(pi/180.0*(Event->alpha+90)) + ycroin*sin(pi/180.0*(Event->alpha+90)));
  //Event->u0 = u0croin*rcroin - (-xcroin*sin(pi/180.0*(Event->alpha+90)) + ycroin*cos(pi/180.0*(Event->alpha+90)));

  //cout << "(t0croin,u0croin)=( " << t0croin << " , " << u0croin << " )\n";
  //cout << "(t0,u0)=( " << Event->t0 << " , " << Event->u0 << " )\n";

  Event->u0max = u0max*rcroin;
  Event->w = Event->u0max*Event->raww / (Event->tE_r/Event->tE_h);
 
}
