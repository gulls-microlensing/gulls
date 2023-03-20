/*! \file 
\brief Functions to create the microlensing event

This file contains functions that computes or collects the parameters required for the microlensing event. 
 */
#include<string>

#include "buildEvent.h"
#include "strfns.h"
//#include "getSynthGalvals.h"
#include "random.h"
#include "structures.h"
#include "constdefs.h"
#include "astroFns.h"
#include "integerPowers.h"

#define DEBUGVAR 0

void addstars(struct event *Event, struct obsfilekeywords World[], 
	      vector<vector<vector<vector<double> > > >* sf, 	      
	      vector<vector<double> >* sfdata, 
	      struct filekeywords *Paramfile, struct slcat* Lenses, 
	      long* idum);
void computeBlending(struct event *Event, struct obsfilekeywords World[],
		     struct filekeywords *Paramfile, struct slcat* Sources);
void drawsl(struct filekeywords* Paramfile, struct event *Event, 
	    struct slcat *Sources, struct slcat *Lenses, long* idum);
void getPlanetvals(struct event* Event, vector<struct pcat>* Planets, struct slcat* Lenses, int sdx);
void compute_u0(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event* Event, long* idum);

void buildEvent(struct event *Event, struct obsfilekeywords World[], 
		vector<vector<vector<vector<double> > > >* starfield, 
		vector<vector<double> >* starfielddata, 
		struct filekeywords *Paramfile, struct slcat *Sources, 
		struct slcat *Lenses, vector<struct pcat> *Planets, 
		int sdx, char* instance_, long *idum)
{
  Event->instance = atoi(instance_);
  Event->id = sdx;

  drawsl(Paramfile, Event, Sources, Lenses, idum);
  getPlanetvals(Event, Planets, Lenses, sdx);

  addstars(Event, World, starfield, starfielddata, Paramfile, Lenses, idum);
  computeBlending(Event, World, Paramfile, Sources);

  compute_u0(Paramfile, World, Event, idum);

  //Place some conditions to speed things up...
      
  double alp = (Event->alpha>=180?360-Event->alpha:Event->alpha);

  if(Event->u0>0.02)
    {
      if(Event->params[SS]<0.1) Event->allsat=1;
      if(Event->params[SS]>7.0 && Event->alpha>100 && Event->alpha<250 && Event->u0>0.5) Event->allsat=1;
      if(log10(Event->params[QQ])<-4.5)
	{
	  if(alp>90.5&&(Event->params[SS]>1.3||Event->params[SS]<0.24)) Event->allsat=1;
	  if(alp<89.5&&Event->params[SS]<0.75) Event->allsat=1;
	  if(Event->alpha<89.5 && Event->params[SS] > 3.05 + 1.0/((90-alp)/55.0)) Event->allsat=1;
	}
    }
  //if(DEBUGVAR) printf("generateParams\n");
  //generateParams(Event,World,Paramfile,idum);
}

void addstars(struct event *Event, struct obsfilekeywords World[], 
	      vector<vector<vector<vector<double> > > >* sf, 	      
	      vector<vector<double> >* sfdata, 
	      struct filekeywords *Paramfile, struct slcat* Lenses, long* idum)
{
  int obsidx, ldx, rdx, sdx;

  //Setup the starfield for the image

  int field = Event->field;
  int filter;
  double solid_angle;
  int xmin,xmax,ymin,ymax,nrepeats;
  int xsub,ysub;
  double Afield;
  double x=0,y=0;
  double pxscl0;
  double Nsub0;
  double convfac;

  pxscl0 = World[0].im.psf.pixscale;
  Nsub0 = double(World[0].im.psf.Nsub);

  //initialize the images
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //reset the detectors
      World[obsidx].im.reset_image();
      World[obsidx].im.reset_detector();

      World[obsidx].im.set_largepsfmag(Paramfile->large_psf_mag);

      int filter = World[obsidx].filter;

      //calculate the position of the source in each image
      if(obsidx==0)
	{
	  //place the star somewhere in the central pixel of the first detector
	  x = (floor(World[obsidx].im.Xpix/2.0) + ran2(idum))*pxscl0;
	  y = (floor(World[obsidx].im.Ypix/2.0) + ran2(idum))*pxscl0;
	}

      //now work out the pixel and subpixel coordinates
      Event->xsub[obsidx] = int(floor(x * World[obsidx].im.psf.Nsub 
				      / World[obsidx].im.psf.pixscale));
      Event->ysub[obsidx] = int(floor(y * World[obsidx].im.psf.Nsub 
				      / World[obsidx].im.psf.pixscale));
      Event->xpix[obsidx] = int(floor(x / World[obsidx].im.psf.pixscale));
      Event->ypix[obsidx] = int(floor(y / World[obsidx].im.psf.pixscale));

      //add the lens star to the image
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       Lenses->mags[Event->lens][filter]);
    }

  //for each field level
  for(ldx=0;ldx<NUM_STARFIELD_LEVELS;ldx++)
    {
      //calculate the required field dimensions
      solid_angle = (*sfdata)[field][ldx];
      //cout << "solidangle = " << solid_angle << endl;
      World[0].im.field_dimensions(solid_angle, xmin, xmax,
				   ymin, ymax, Afield, nrepeats);
      
      int nstars = int((*sf)[field][ldx].size());
      //cout << "nstars = " << nstars << endl;
      //cout << "nrepeats = " << nrepeats << endl;

      //for each repetition
      for(rdx=0;rdx<nrepeats;rdx++)
	{
	  //for each star
	  for(sdx=0;sdx<nstars;sdx++)
	    {
	      //choose a random position
	      x = (xmin+(xmax-xmin)*ran2(idum)) * pxscl0/Nsub0;
	      y = (ymin+(ymax-ymin)*ran2(idum))	* pxscl0/Nsub0;
 
	      //for each observatory
	      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  convfac = double(World[obsidx].im.psf.Nsub) 
		    / World[obsidx].im.psf.pixscale;
		  filter = World[obsidx].filter;

		  //work out where the star is in subpixel coordinates  
		  xsub = int(floor(x*convfac));
		  ysub = int(floor(y*convfac));
		    
		  World[obsidx].im.addstar(xsub, ysub, 
					   (*sf)[field][ldx][sdx][filter]);

		  if(DEBUGVAR && (*sf)[field][ldx][sdx][filter]<13 && xsub>0 && xsub<World[obsidx].im.Xpix*World[obsidx].im.psf.Nsub && ysub>=0 && ysub<World[obsidx].im.Ypix*World[obsidx].im.psf.Nsub)
		    {
		      cerr << "FLAG mag " << (*sf)[field][ldx][sdx][filter] << " star at pixel " << xsub/World[obsidx].im.psf.Nsub+1 << "," << ysub/World[obsidx].im.psf.Nsub+1 << " in band " << obsidx << endl;
		    }
		} //end for each observatory

	    } //end for each star
	  
	} //end for each repetition

      //add the background due to PSF tails
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	{
	  double scalefac = pxscl0*World[obsidx].im.psf.Nsub 
	    / (World[obsidx].im.psf.pixscale*Nsub0);
	  World[obsidx].im.addpsfbg(int(scalefac*xmin), int(scalefac*xmax), 
				    int(scalefac*ymin), int(scalefac*ymax));
	}

    } //end for each field level
      
}

void computeBlending(struct event *Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources)
{
  int filter;
  double nc, err, nc0;
  int satflag;
  int obsidx;
  int sn = Event->source;

  Event->allsat=1;

  //Work out the blending and if the event is always saturated
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    { 
      filter = World[obsidx].filter;

      Event->allsatobs[obsidx]=1;  //initialize

      //add the image background
      //Actually no, calculate blending in the no background case
      //World[obsidx].im.addbg();

      //blending fraction: calculate the flux without the source
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx],
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&nc, &err, &satflag, false);

      //add the baseline source and recalculate blending
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       Sources->mags[sn][filter]);
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&nc0, &err, &satflag, false);
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       Sources->mags[sn][filter]);

      
      if(!satflag) 
	{
	  Event->allsat=0;
	  Event->allsatobs[obsidx]=0;
	}
	      
      //work out the blending fraction
      Event->fs[obsidx] = (nc0-nc)/nc0;
      Event->baselineFlux[obsidx] = nc0/World[obsidx].mintexp;

      Event->peakpoint=0;

      //subtract the background if neccessary
      //if(!CONSTANT_BACKGROUND) World[obsidx].im.subbg();
    }
}

void drawsl(struct filekeywords* Paramfile, struct event *Event, struct slcat *Sources, struct slcat *Lenses, long* idum)
{
  //draw a lens and a source from their catalogues, make sure its a valid pair 
  //and compute the relevant quantities

  int sn, ln, fn; //source, lens and field numbers

  //choose a random field
  fn = Paramfile->validFields[randint(0,Paramfile->validFields.size()-1, 
				      idum)];
  
  do
    {
      //choose a source and lens
      sn = randint(Sources->start[fn], Sources->end[fn], idum);
      ln = randint(Lenses->start[fn], Lenses->end[fn], idum);

    } while(Sources->data[sn][DIST] <= Lenses->data[ln][DIST]
	    || (Sources->data[sn][MUL] == Lenses->data[ln][MUL] 
		&& Sources->data[sn][MUB] == Lenses->data[ln][MUB]));

  //store the choice
  Event->field = fn;
  Event->source = sn;
  Event->lens = ln;

  //calculate its fundamental microlensing properties
    
  //double tE, rE, thE, piE, rs, mu, vt;

  //relative ls proper motion - lens primary is the origin
  //may need some cos(b)'s in here somewhere???
  //in mas/yr

  //einstein radius
  double x = Lenses->data[ln][DIST]/Sources->data[sn][DIST];
  //cout << x << endl;

  //in AU
  Event->rE = rEsun * sqrt(Lenses->data[ln][MASS] 
			   * Sources->data[sn][DIST] * (1-x) * x);

  //in mas
  Event->thE = Event->rE/Lenses->data[ln][DIST];

  Event->piE = (1-x)/Event->rE;

  /*double ul,vl,wl,us,vs,ws,ll,lb,sl,sb,dl,ds;

  ul = (Lenses->data[ln][UU]);       us = (Sources->data[sn][UU]);
  vl = (Lenses->data[ln][VV]);       vs = (Sources->data[sn][VV]);
  wl = (Lenses->data[ln][WW]);       ws = (Sources->data[sn][WW]);
  ll = Lenses->data[ln][LL]*TO_RAD;  lb = Lenses->data[ln][BB]*TO_RAD;
  sl = Sources->data[sn][LL]*TO_RAD; sb = Sources->data[sn][BB]*TO_RAD;
  dl = Lenses->data[ln][DIST];       ds = Sources->data[sn][DIST];

  double vll = ul*sin(ll)*cos(lb) + vl*cos(ll)*cos(lb) + wl*sin(lb);
  double vlb = ul*sin(lb) + vl*sin(ll)*cos(lb) + wl*cos(lb);
  double vsl = us*sin(sl)*cos(sb) + vs*cos(sl)*cos(sb) + ws*sin(sb);
  double vsb = us*sin(sb) + vs*sin(sl)*cos(sb) + ws*cos(sb);*/

  //in mas yr-1
  Event->murel = sqrt(sqr(Sources->data[sn][MUL]-Lenses->data[ln][MUL]) + sqr(Sources->data[sn][MUB]-Lenses->data[ln][MUB]));

  //in kms-1
  Event->vt = Event->murel * Lenses->data[ln][DIST] * AU/1000 / SECINYR;

  //in km s-1
  //Event->vt = sqrt(sqr((vsl*dl/ds-vll)*cos(lb)) + sqr(vsb*dl/ds-vlb));

  //in mas yr-1
  //Event->murel = sqrt(sqr(vsl/ds - vll/dl) + sqr(vsb/ds - vlb/dl))
  //  * 1000*SECINYR/AU;

  //in days
  Event->tE = DAYINYR * Event->thE / Event->murel;

  //in Einstein radii
  Event->rs = (Sources->data[sn][RADIUS] * Rsun / Sources->data[sn][DIST]) / Event->thE;
  //radius (Rsun) -> AU / Ds (kpc) -> mas / thetaE (mas) = ratio

  //rate weighting
  Event->raww = Event->rE * Event->vt;

  //positions
  //Event->l = Lenses->data[ln][LL];
  //Event->b = Lenses->data[ln][BB];
  Event->l = Lenses->l[fn] + (ran2(idum)-0.5)*Lenses->dl[fn];
  Event->b = Lenses->b[fn] + (ran2(idum)-0.5)*Lenses->db[fn];
  eq2gal(Event->l*TO_RAD, Event->l*TO_RAD, 'g', &Event->ra, &Event->dec);
  if(Event->ra<0) Event->ra += 2*PI;
  //may want to convert ra and dec to degrees here?

  //The random parameters

  Event->t0 = double(NUM_SIM_DAYS)*ran2(idum);
  Event->alpha = 360.0 * ran2(idum);

  //u0 will be calculated after we know the blending
}

void getPlanetvals(struct event* Event, vector<struct pcat>* Planets, struct slcat* Lenses, int sdx)
{
  //extract and calculate the planet parameters

  int ln = Event->lens;
  Event->params.resize(NPLANETINPUT+NPLANETDERIV);

  //extract the input data
  for(int i=0;i<NPLANETINPUT;i++)
    {
      Event->params[i] = (*Planets)[sdx].data[i];
    }

  //Semimajor axis is currently holding the period - swap them...
  Event->params[TT] = Event->params[AA];
  Event->params[AA] = pow(
			  (Event->params[PMASS] + Lenses->data[ln][MASS]) 
			  * sqr(Event->params[TT]),1.0/3.0);

  //Calculate the derived planet properties
  Event->params[QQ] = Event->params[PMASS] / Lenses->data[ln][MASS];

  double sqrt1pq = sqrt(1+Event->params[QQ]);

  Event->rE *= sqrt1pq;
  Event->tE *= sqrt1pq;
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

}

void compute_u0(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event* Event, long* idum)
{

  double umaxmax=3.0;
  double umaxmin=0.01;

  int obsidx = Paramfile->principle_observatory;

  //Calculating u0

  //We require a change in flux of 2% @ peak to detect it
  double Amin = Paramfile->Amin;

  //Compute using the total blending in primary band (lens+blend)

  double fs = Event->fs[obsidx];
  double mumin = (Amin-1+fs)/fs;

  //ensure that the microlensing event can be seen, but don't go too far
  //into the finite source regime
  double u0max = sqrt(2.0*sqrt(1.0+1.0/(mumin*mumin-1.0)) - 2.0);
  u0max = (u0max>umaxmax?umaxmax:u0max);
  u0max = (u0max<umaxmin?umaxmin:u0max);

  Event->u0 = u0max*ran2(idum);

  Event->w = u0max*Event->raww;
  Event->u0max = u0max;
}
