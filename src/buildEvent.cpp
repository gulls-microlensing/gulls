/*! \file 
\brief Functions to create the microlensing event

This file contains functions that computes or collects the parameters required for the microlensing event. 
 */
#include<string>

#include "buildEvent.h"
#include "strfns.h"
#include "random.h"
#include "structures.h"
#include "constdefs.h"
#include "astroFns.h"
#include "ephem.h"
#include "getPlanetvals.h"
#include "photometryCodes.h"

#define DEBUGVAR 0

void buildEvent(struct event *Event, struct obsfilekeywords World[], 
		vector<vector<vector<double> > >* starfield, 
		vector<double>* starfielddata, 
		struct filekeywords *Paramfile, struct slcat *Sources, 
		struct slcat *Lenses, int sdx, char* instance_, long *idum)
{
  Event->instance = atoi(instance_);
  Event->id = sdx;

  //clear the data vector
  Event->data.clear();

  //Set up obsgroups
  if(int(Event->obsgroups.size())==0)
    {
      if(Paramfile->verbosity>1) cout << "drawsl" << endl;
      setupObsGroups(Paramfile, Event);
    }

  if(Paramfile->verbosity>1) cout << "drawsl" << endl;
  drawsl(Paramfile, World, Event, Sources, Lenses, idum);
  //getPlanetvals(Event, Planets, Lenses, sdx);

  if(Paramfile->verbosity>1) cout << "addstars" << endl;
  addstars(Event, World, starfield, starfielddata, Paramfile, Sources, Lenses, idum);
  if(Paramfile->verbosity>1) cout << "computeBlending" << endl;
  computeBlending(Event, World, Paramfile, Sources, Lenses);

  if(Paramfile->verbosity>1) cout << "compute_u0" << endl;
  compute_u0(Paramfile, World, Event, idum);

  Event->peakpoint=0;

  if(ran2(Paramfile->seed)<Paramfile->outputLightcurve) Event->outputthis=1;
  else Event->outputthis=0;
  cout << "outputLightcurve = " << Paramfile->outputLightcurve << " " << Event->outputthis << endl;
}

void addstars(struct event *Event, struct obsfilekeywords World[], 
	      vector<vector<vector<double> > >* sf, 	      
	      vector<double>* sfdata, 
	      struct filekeywords *Paramfile, struct slcat *Sources, 
	      struct slcat* Lenses, long* idum)
{
  int obsidx, ldx, rdx, sdx;

  //Setup the starfield for the image

  //int field = Event->field;
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
  Event->sl.clear();
  Event->sl.resize(Paramfile->numobservatories);

  //initialize the images
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //reset the detectors
      World[obsidx].im.reset_image();
      World[obsidx].im.reset_detector();
      World[obsidx].ref.reset_image();
      World[obsidx].ref.reset_detector();

      //reset the starlists
      Event->sl[obsidx].reset();

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
      if(Paramfile->lenslight)
	{
	  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   Lenses->mags[Event->lens][filter]);
	}
    }

  //for each field level
  for(ldx=0;ldx<NUM_STARFIELD_LEVELS;ldx++)
    {
      //calculate the required field dimensions
      solid_angle = (*sfdata)[ldx];
      //World[0].im.field_dimensions(solid_angle, xmin, xmax,
      //				   ymin, ymax, Afield, nrepeats);
      World[0].im.field_dimensions(solid_angle, xmin, xmax,
      				   ymin, ymax, Afield);
      
      int nstars = int((*sf)[ldx].size());
      if(nstars==0) break;
      int npick = poisson(nstars*Afield/solid_angle,idum);

      if(Paramfile->verbosity>2)
	cout << "nstars, npick, solidangle, xmin, xmax, ymin, ymax, (xmax-xmin), (ymax-ymin), Afield, mean " << nstars << " " << npick << " " << solid_angle << " " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << (xmax-xmin) << " " << (ymax-ymin) << " " << Afield << " " << nstars*Afield/solid_angle << endl;

      //for each repetition
      //for(rdx=0;rdx<nrepeats;rdx++)
      //{
      //for each star
      //for(sdx=0;sdx<nstars;sdx++)
      for(int star=0;star<npick;star++)
	{
	  //pick a random star
	  sdx = randint(0,nstars-1,idum);
	  
	  //choose a random position
	  x = (xmin+(xmax-xmin)*ran2(idum)) * pxscl0/Nsub0;
	  y = (ymin+(ymax-ymin)*ran2(idum)) * pxscl0/Nsub0;
 
	  //for each observatory
	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      convfac = double(World[obsidx].im.psf.Nsub) 
		/ World[obsidx].im.psf.pixscale;
	      filter = World[obsidx].filter;

	      //work out where the star is in subpixel coordinates  
	      xsub = int(round(x*convfac));
	      ysub = int(round(y*convfac));

	      //cout << xsub << " " << ysub << " " << convfac << endl;
		    
	      World[obsidx].im.addstar(xsub, ysub, 
				       (*sf)[ldx][sdx][filter],
				       &Event->sl[obsidx]);

	      if(Paramfile->verbosity>1 && (*sf)[ldx][sdx][filter]<13 && xsub>0 && xsub<World[obsidx].im.Xpix*World[obsidx].im.psf.Nsub && ysub>=0 && ysub<World[obsidx].im.Ypix*World[obsidx].im.psf.Nsub)
		{
		  cerr << "FLAG mag " << (*sf)[ldx][sdx][filter] << " star at pixel " << xsub/World[obsidx].im.psf.Nsub+1 << "," << ysub/World[obsidx].im.psf.Nsub+1 << " in band " << obsidx << endl;
		}
	    } //end for each observatory

	} //end for each star
	  
      //} //end for each repetition

    } //end for each field level
      
}

void computeBlending(struct event *Event, struct obsfilekeywords World[], struct filekeywords *Paramfile, struct slcat* Sources, struct slcat* Lenses)
{
  int filter;
  int satflag=0;
  int obsidx;
  int sn = Event->source;
  int tmpsatflag=0;
  //int ln = Event->lens;

  vector<double> phot0, phot1; //photometry before and after adding source
  
  double allbg; //the magnitude of all combined backgrounds

  Event->allsat=1;

  //Work out the blending and if the event is always saturated
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    { 
      filter = World[obsidx].filter;

      Event->allsatobs[obsidx]=1;  //initialize

      //add the image background
      if(Paramfile->verbosity>2) 
	cout << __FILE__ << " " << __FUNCTION__ << ": Add background " << World[obsidx].zodiflux[0]
	     << endl;
      allbg = 20.0 - 2.5*log10(World[obsidx].constbackground 
			       + World[obsidx].zodiflux[0]
			       + pow(10,-0.4*(World[obsidx].skybackground-20))
			       );
      World[obsidx].im.set_background(allbg);
      World[obsidx].im.addbg();

      //blending fraction: calculate the flux without the source
      if(Paramfile->verbosity>2) 
	cout << __FILE__ << " " << __FUNCTION__ << ": Photometry w/o source"
	     << endl;
      World[obsidx].im.wis_photometry(Event->xsub[obsidx], 
				      Event->ysub[obsidx],
				      World[obsidx].mintexp, 1,
				      &phot0, &tmpsatflag);
      satflag |= tmpsatflag;

      //add the baseline source and recalculate blending
      if(Paramfile->verbosity>2) 
	cout <<  __FILE__ << " " << __FUNCTION__ << ": Add source"
	     << endl;
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       Sources->mags[sn][filter]);

      if(Paramfile->verbosity>2) 
	cout << __FILE__ << " " << __FUNCTION__ << ": Photometry with source"
	     << endl;
      World[obsidx].im.wis_photometry(Event->xsub[obsidx], 
				      Event->ysub[obsidx],
				      World[obsidx].mintexp, 1,
				      &phot1, &tmpsatflag);

      //Store the blending fraction
      if(World[obsidx].photcode<2) //aperture photometry
	{
	  Event->fs[obsidx] = (phot1[0]-phot0[0])/phot1[0];
	  Event->baselineFlux[obsidx] = phot1[0]/World[obsidx].mintexp;
	}
      else //weighted photometry
	{
	  Event->fs[obsidx] = (phot1[4]-phot0[4])/phot1[4];
	  Event->baselineFlux[obsidx] = phot1[4]/World[obsidx].mintexp;
	}


      //Check for saturation due to bleeding by generating an exposure

      double dummycounts, dummyerr;
      World[obsidx].im.reset_detector();
      World[obsidx].im.expose(World[obsidx].mintexp); //assume all exposures same length
      World[obsidx].im.photometry(Event->xpix[obsidx], 
				  Event->ypix[obsidx],&dummycounts,&dummyerr,&tmpsatflag);
      World[obsidx].im.reset_detector();

      if(satflag==0&&tmpsatflag!=0&&Paramfile->verbosity>0)
	{
	  cout << "All sat flag set due to bleeding test" << endl;
	}

      satflag |= tmpsatflag;

      if(!satflag) 
	{
	  Event->allsat=0;
	  Event->allsatobs[obsidx]=0;
	}    

      //remember to remove the background and the source again
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       Sources->mags[sn][filter]);

      if(World[obsidx].photcode<0) //setup the fast photometry
	{
	  World[obsidx].im.setup_fast_photometry(Event->xsub[obsidx], 
						 Event->ysub[obsidx], 
						 Event->xpix[obsidx], 
						 Event->ypix[obsidx], 
						 Sources->mags[sn][filter], 
						 World[obsidx].exptime[0],
						 World[obsidx].nstack[0]);

	  //fast_blend includes all the counts from the unmagnified source too - well named past self!
	  Event->fs[obsidx] = World[obsidx].im.fast_src/(World[obsidx].im.fast_blend);
	  Event->baselineFlux[obsidx] = (World[obsidx].im.fast_blend)/(World[obsidx].exptime[0]*World[obsidx].nstack[0]);
	}

      World[obsidx].im.subbg();

    }
}

void drawsl(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event *Event, struct slcat *Sources, struct slcat *Lenses, long* idum)
{
  //draw a lens and a source from their catalogues, make sure its a valid pair 
  //and compute the relevant quantities

  int sn, ln; //source, lens and field numbers
  vector<double> lb(2); //galactic coordinates of the event
  double x;
  vector<double> pmgal(2);


  //choose a random field - now obsolete
  //fn = Paramfile->validFields[randint(0,Paramfile->validFields.size()-1, 
  //			      idum)];
  
  do
    {
      //choose a source and lens
      sn = randint(Sources->start, Sources->end, idum);
      ln = randint(Lenses->start, Lenses->end, idum);

    } while(Sources->data[sn][DIST] <= Lenses->data[ln][DIST]
	    || (Sources->data[sn][MUL] == Lenses->data[ln][MUL] 
		&& Sources->data[sn][MUB] == Lenses->data[ln][MUB]));

  //store the choice
  Event->field = Paramfile->choosefield;
  Event->source = sn;
  Event->lens = ln;

  //fractional lens source distance
  x = Lenses->data[ln][DIST]/Sources->data[sn][DIST];

  //positions
  //randomly choose an l,b somewhere in the box
  Event->l = Lenses->l + (ran2(idum)-0.5)*Lenses->dl;
  Event->b = Lenses->b + (ran2(idum)-0.5)*Lenses->db;
  lb[0] = Event->l*TO_RAD;
  lb[1] = Event->b*TO_RAD;
  eq2gal(lb[0], lb[1], 'g', &Event->ra, &Event->dec);
  if(Event->ra<0) Event->ra += 2*PI;

  //The random parameters

  Event->t0 = double(Paramfile->NUM_SIM_DAYS)*ran2(idum);
  Event->alpha = 360.0 * ran2(idum);

  //u0 will be calculated after we know the blending

  //calculate the fundamental microlensing properties
    
  //double tE, rE, thE, piE, rs, mu, vt;

  //einstein radius
  //cout << x << endl;

  //in AU
  Event->rE = rEsun * sqrt(Lenses->data[ln][MASS] 
			   * Sources->data[sn][DIST] * (1-x) * x);

  //in mas
  Event->thE = Event->rE/Lenses->data[ln][DIST];

  //relative ls proper motion - lens motion relative to the source
  //in mas/yr
  //calculate the heliocentric relative proper motion
  pmgal[0] = Lenses->data[ln][MUL]-Sources->data[sn][MUL];
  pmgal[1] = Lenses->data[ln][MUB]-Sources->data[sn][MUB];

  //work out its absolute value
  Event->murel = qAdd(pmgal[0],pmgal[1]);

  //calculate the parallax
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

  //in kms-1
  Event->vt = Event->murel * Lenses->data[ln][DIST] * AU/1000 / SECINYR;

  //in km s-1
  //Event->vt = sqrt(sqr((vsl*dl/ds-vll)*cos(lb)) + sqr(vsb*dl/ds-vlb));

  //in mas yr-1
  //Event->murel = sqrt(sqr(vsl/ds - vll/dl) + sqr(vsb/ds - vlb/dl))
  //  * 1000*SECINYR/AU;

  //in days
  Event->tE_h = DAYINYR * Event->thE / Event->murel;

  //in Einstein radii
  Event->rs = (Sources->data[sn][RADIUS] * Rsun / Sources->data[sn][DIST]) / Event->thE;
  //radius (Rsun) -> AU / Ds (kpc) -> mas / thetaE (mas) = ratio

  //rate weighting
  Event->raww = Event->thE * Event->murel;

}

void compute_u0(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event* Event, long* idum)
{

  /*
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
  */

  Event->u0max = Paramfile->u0max;
  Event->u0 = Event->u0max*(2*ran2(idum)-1);
  Event->w = Event->u0max*Event->raww;
  
}

//not actually used here, but provided as a utility. May include it as a 
//parameter file option later, and build in the coping mechanisms in the 
//analysis scripts
int inSeason(double t0, struct filekeywords* Paramfile, struct obsfilekeywords World[])
{
  if(int(floor(t0))<0 || int(floor(t0))>Paramfile->NUM_SIM_DAYS) return 0;
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    { 
      if(World[obsidx].weatherSequence[int(floor(t0*4))]>0) return 1;
    }
  return 0;
}

void setupParallax(double tref, struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event *Event, struct slcat *Sources, struct slcat *Lenses)
{
  int sn = Event->source;
  int ln = Event->lens;

  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //Prepare the epochs correctly
      vector<double> jdepochs(1,Paramfile->simulation_zerotime);

      if(int(jdepochs.size())>0)
	{
	  Event->pllx[obsidx].reset();
	  Event->pllx[obsidx].set_reference(Paramfile->simulation_zerotime+tref,&World[0].orbit);
	  Event->pllx[obsidx].set_orbit(&World[obsidx].orbit);
	  Event->pllx[obsidx].set_lb(Event->l, Event->b);
	  Event->pllx[obsidx].set_pm_lb(Lenses->data[ln][MUL]-Sources->data[sn][MUL],Lenses->data[ln][MUB]-Sources->data[sn][MUB]);
	  Event->pllx[obsidx].set_piE(Event->piE);
	  Event->pllx[obsidx].set_tE_h(Event->tE_h);

	  //Do this for a dummy epoch now, but redo this at the end of 
	  //timeSequencer
	  Event->pllx[obsidx].load_epochs(&World[obsidx].jd);
	  Event->pllx[obsidx].initialize();

	  if(Paramfile->verbosity>3)
	    {
	      for(int idx=0;idx<int(World[obsidx].jd.size());idx++)
		{
		  printf("NEShift event %d obs %d %14.6f %g %g\n",Event->id,obsidx,World[obsidx].jd[idx],Event->pllx[obsidx].NEshift[idx][0] + (Event->pllx[obsidx].epochs[idx]-Event->pllx[obsidx].tref)*Event->pllx[obsidx].vref[0],Event->pllx[obsidx].NEshift[idx][1] + (Event->pllx[obsidx].epochs[idx]-Event->pllx[obsidx].tref)*Event->pllx[obsidx].vref[1]);
		}
	    }
	}

      if(Paramfile->verbosity>0)
	{
	  cerr << "Observatory " << obsidx << ":" << endl;
	  cerr << "Orbit: " << endl;
	  cerr << World[obsidx].orbit[0].xh << " " << World[obsidx].orbit[0].yh << " " << World[obsidx].orbit[0].zh << endl;
	    
	  cerr << "tref = " << Event->pllx[obsidx].tref << endl;
	  cerr << "l,b,a,d = " << Event->pllx[obsidx].l << " " << Event->pllx[obsidx].b << " " << Event->pllx[obsidx].a << " " << Event->pllx[obsidx].d << " " << endl;
	  cerr << "mul, mub, mua, mud = " << Event->pllx[obsidx].mul << " " << Event->pllx[obsidx].mub << " " << Event->pllx[obsidx].mua << " " << Event->pllx[obsidx].mud << " " << endl;
	  cerr << "piEN, piEE, piEll, piErp, piE = " << Event->pllx[obsidx].piEN << " " << Event->pllx[obsidx].piEE << " " << Event->pllx[obsidx].piEll << " " << Event->pllx[obsidx].piErp << " " << Event->pllx[obsidx].piE << endl;
	  cerr << "tE_h, tE_r = " << Event->pllx[obsidx].tE_h << " " << Event->pllx[obsidx].tE_r << endl;
	}

    }
  
  Event->piEN = Event->pllx[0].piEN*Event->piE;
  Event->piEE = Event->pllx[0].piEE*Event->piE;
  Event->tE_h = Event->pllx[0].tE_h;
  Event->tE_r = Event->pllx[0].tE_r;
}

void setupObsGroups(struct filekeywords *Paramfile, struct event *Event)
{
  
  //Parse the observatory groups string
  string::size_type start=0;
  string::size_type end;
  //string rest=string(Paramfile->obsgroupstr);
  //string("(ALL)(2,3,1)PERMUTE_PAIRS(ALL,1,2)PERMUTE_REMOVE,PERMUTE_PAIRS(0,1,ALL)");
  string rest=string(Paramfile->obsgroupstr);
  string thisgrp;
  vector<int> tmp;
  int obsidx;

  while(start!=string::npos)
    {
      //groups are contained within parentheses, unless special codes
      start = rest.find_first_of("("); 
      //should be 0 unless a special multi-group codeword
      if(Paramfile->verbosity>1) 
	cout << "setupObsGroups: start = " << start << endl;

      if(start>0||start==string::npos)
	{
	  string special=rest.substr(0,start);
	  if(special.find("PERMUTE_PAIRS")!=string::npos)
	    {
	      vector<int> pair(2);
	      for(obsidx=0; obsidx<Paramfile->numobservatories; obsidx++)
		{
		  for(int obsjdx=obsidx+1; obsjdx<Paramfile->numobservatories; obsjdx++)
		    {
		      pair[0]=obsidx; pair[1]=obsjdx;
		      Event->obsgroups.push_back(pair);
		    }
		}
	    }
	  if(special.find("PERMUTE_REMOVE")!=string::npos)
	    {
	      for(obsidx=Paramfile->numobservatories-1;obsidx>=0;obsidx--)
		{
		  vector<int> set;
		  for(int obsjdx=0;obsjdx<Paramfile->numobservatories;obsjdx++)
		    {
		      if(obsjdx!=obsidx) set.push_back(obsjdx);
		    }
		  Event->obsgroups.push_back(set);
		}
	    }
	  if(special.find("EACH_INDIVIDUAL")!=string::npos)
	    {
	      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  vector<int> set(1,obsidx);
		  Event->obsgroups.push_back(set);
		}
	    }
	}
      if(start!=string::npos)
	{
	  rest = rest.substr(start+1);
	  end = rest.find_first_of(")");
	  thisgrp = rest.substr(0,end-start);
	  Event->obsgroups.push_back(tmp);

	  if(Paramfile->verbosity>1) 
	    {
	      cout << "setupObsGroups: start, rest, end, thisgrp" << endl;
	      cout << "setupObsGroups: " << start << "\t" << rest << "\t" << end << "\t" << thisgrp << endl;
	    }

	  //elements are separated by commas
	  string::size_type elstart=-1;
	  string::size_type comma;
	  do
	    {
	      comma=thisgrp.find_first_of(",",elstart+1);
	      string elstr;
	      if(comma==string::npos) elstr=thisgrp.substr(elstart+1,comma);
	      else elstr = thisgrp.substr(elstart+1,comma-elstart-1);

	      if(Paramfile->verbosity>1) 
		{
		  cout << "setupObsGroups: elstart, comma, elstr" << endl;
		  cout << "setupObsGroups: " << elstart << "\t" << comma << "\t" << elstr << endl;
		}

	      elstart=comma;
	      if(elstr.find("ALL")!=string::npos)
		{
		  Event->obsgroups.back().clear();
		  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		    {
		      Event->obsgroups.back().push_back(obsidx);
		    }
		  elstart=string::npos;
		}
	      else
		{
		  int candidate=atoi(elstr.c_str());
		  if((candidate>0 && candidate<Paramfile->numobservatories) || 
		     (candidate==0 && elstr.find("0")!=string::npos))
		    {
		      Event->obsgroups.back().push_back(candidate);
		    }
		  else
		    {
		      cerr << "Invalid observatory number in group " <<Event->obsgroups.size()-1 << " (" << elstr << ")" << endl;
		      exit(1);
		    }
		}
	    }
	  while(elstart!=string::npos);

	  rest = rest.substr(end+1);
	}
    }

  if(Paramfile->verbosity) 
    {
      cout << "Observing groups specified:" << endl;
      for(int i=0;i<int(Event->obsgroups.size());i++)
	{
	  cout << "\tGroup " << i << endl << "\t\t";
	  for(int j=0;j<int(Event->obsgroups[i].size());j++)
	    {
	      cout << Event->obsgroups[i][j] << " ";
	    }
	  cout << endl; 
	}
    }

  //Setup the size of anything that depends on the groups
  Event->PSPL.resize(Event->obsgroups.size());
  Event->FSPL.resize(Event->obsgroups.size());
  Event->flatchi2.resize(Event->obsgroups.size());
  Event->flag_needFS.resize(Event->obsgroups.size());
  Event->flatlc.resize(Event->obsgroups.size());
  Event->obsgroupoutput.clear();
  Event->obsgroupoutput.resize(Event->obsgroups.size(),string(""));
}

