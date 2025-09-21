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
#include "coords.h"

#define DEBUGVAR 0

void buildEvent(struct event *Event, struct obsfilekeywords World[], 
		vector<vector<vector<double> > >* starfield, 
		vector<double>* starfielddata, 
		struct filekeywords *Paramfile, struct slcat *Sources, 
		struct slcat *Lenses, int sdx, string instance, long *idum)
{
  Event->instance = atoi(instance.c_str());
  Event->id = sdx;

  Event->gamma = Paramfile->LD_GAMMA;

  //clear the data vectors
  Event->data.clear();
  Event->scomp_rs.clear();
  Event->scomp_s.clear();
  Event->scomp_alpha.clear();
  Event->scomp_inc.clear();
  Event->scomp_phase.clear();
  Event->scomp_fsofs1.clear();

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

      //add the lens star and lens and source companions to the image
      if(Paramfile->lenslight)
	{
	  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   Lenses->mags[Event->lens][filter]);
	  if(Paramfile->multiple_lenses)
	    {
	      for(auto lc : Event->lcompanions)
		{
		  if(Paramfile->verbosity>2) cout << "Adding lens companion star (mag=" << Lenses->mags[lc][filter] << ")" << endl; 
		  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
					   Lenses->mags[lc][filter]);
		}
	    }	  
	}
      if(Paramfile->multiple_sources)
	{
	  for(auto sc : Event->scompanions)
	    {
	      if(Paramfile->verbosity>2) cout << "Adding source companion star (mag=" << Sources->mags[sc][filter] << ")" << endl; 
	      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   Sources->mags[sc][filter]);
	    }
	}
    }

  //for each field level
  for(ldx=0;ldx<NUM_STARFIELD_LEVELS;ldx++)
    {
      //calculate the required field dimensions
      solid_angle = (*sfdata)[ldx]; //in square arcsec
      //World[0].im.field_dimensions(solid_angle, xmin, xmax,
      //				   ymin, ymax, Afield, nrepeats);
      World[0].im.field_dimensions(solid_angle, xmin, xmax,
      				   ymin, ymax, Afield);
      
      int nstars = int((*sf)[ldx].size());
      if(nstars==0) continue;
      int npick = poisson(nstars*Afield/solid_angle,idum);

      if(Paramfile->verbosity>2)
	cout << "nstars, npick, solidangle(arcsec^2), xmin, xmax, ymin, ymax, (xmax-xmin), (ymax-ymin), Afield, mean " << nstars << " " << npick << " " << solid_angle << " " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << (xmax-xmin) << " " << (ymax-ymin) << " " << Afield << " " << nstars*Afield/solid_angle << endl;

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
		{
		  cout << __FILE__ << " " << __FUNCTION__ << ": Add background " << World[obsidx].zodiflux[0]
			   << endl;
		}
      allbg = 20.0 - 2.5*log10(World[obsidx].constbackground 
			       + World[obsidx].zodiflux[0]
			       + pow(10,-0.4*(World[obsidx].skybackground-20))
			       );
      World[obsidx].im.set_background(allbg);
      World[obsidx].im.addbg();

      //blending fraction: calculate the flux without the source
      if(Paramfile->verbosity>2)
		{
		  cout << __FILE__ << " " << __FUNCTION__ << ": Photometry w/o source"
			   << endl;
		}
      World[obsidx].im.wis_photometry(Event->xsub[obsidx], 
									  Event->ysub[obsidx],
									  World[obsidx].mintexp, 1,
									  &phot0, &tmpsatflag);
      satflag |= tmpsatflag;

      //add the baseline source and recalculate blending
      if(Paramfile->verbosity>2)
		{
		  cout <<  __FILE__ << " " << __FUNCTION__ << ": Add source"
			   << endl;
		}
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
							   Sources->mags[sn][filter], false, true);
	  
      if(Paramfile->verbosity>2)
		{
		  cout << __FILE__ << " " << __FUNCTION__ << ": Photometry with source"
			   << endl;
		}
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
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
							   Sources->mags[sn][filter], true, true);
	  
      if(World[obsidx].photcode<0) //setup the fast photometry
	{
	  World[obsidx].im.setup_fast_photometry(Event->xsub[obsidx], 
						 Event->ysub[obsidx], 
						 Event->xpix[obsidx], 
						 Event->ypix[obsidx], 
						 Sources->mags[sn][filter], 
						 World[obsidx].mintexp,
						 1.0);
												 //World[obsidx].exptime[0],
												 //World[obsidx].nstack[0]);
		  
	  //fast_blend includes all the counts from the unmagnified source too - well named past self!
	  Event->fs[obsidx] = World[obsidx].im.fast_src/(World[obsidx].im.fast_blend);
	  //Event->baselineFlux[obsidx] = (World[obsidx].im.fast_blend)/(World[obsidx].exptime[0]*World[obsidx].nstack[0]);
	  Event->baselineFlux[obsidx] = (World[obsidx].im.fast_blend)/(World[obsidx].mintexp*1.0);
	  cout << Event->fs[obsidx] << " " << Event->baselineFlux[obsidx] << endl;
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

    } while(Sources->data[sn][Sources->DIST] <= Lenses->data[ln][Lenses->DIST]
	    || (Sources->data[sn][Sources->MUL] == Lenses->data[ln][Lenses->MUL] 
		&& Sources->data[sn][Sources->MUB] == Lenses->data[ln][Lenses->MUB]));

  //store the choice
  Event->field = Paramfile->choosefield;
  Event->source = sn;
  Event->lens = ln;

  //fractional lens source distance
  x = Lenses->data[ln][Lenses->DIST]/Sources->data[sn][Sources->DIST];

  //positions
  //randomly choose an l,b somewhere in the box
  Event->l = Lenses->l + (ran2(idum)-0.5)*Lenses->dl;
  Event->b = Lenses->b + (ran2(idum)-0.5)*Lenses->db;
  lb[0] = Event->l*TO_RAD;
  lb[1] = Event->b*TO_RAD;
  eq2gal(lb[0], lb[1], 'g', &Event->ra, &Event->dec);
  if(Event->ra<0) Event->ra += 2*PI;

  //Handle multiplicity
  Event->scompanions.clear();
  if(Paramfile->multiple_sources>0)
    {
      int isbinary = Sources->data[sn][Sources->datadict["Is_Binary"]];
      if(isbinary>0)
	{
	  if(Paramfile->verbosity>2) cout << "Source " << sn << " is a multiple ";
	  //Currently set up so that companions immediately trail the primary in the catalog
	  if(isbinary==1)
	    {
	      if(Paramfile->verbosity>2) cout << "and is the primary." << endl;
	      for(int i=sn+1;i<Sources->data.size();i++)
		{
		  if(Sources->data[i][Sources->datadict["primary_ID"]]==Sources->data[sn][Sources->datadict["ID"]])
		    {
		      if(Paramfile->verbosity>2) cout << "Adding star " << i << " to source system." << endl;
		      Event->scompanions.push_back(i);
		      //add companion properties to the event data here?
		    }
		  else
		    {
		      break;
		    }
		}
	    } //end isbinary==1
	  else if(isbinary>=2)
	    {
	      if(Paramfile->verbosity>2) cout << "and is a companion." << endl;
	      //find the primary id's position
	      int primarysn=-1;
	      for(int i=0;i<Sources->data.size();i++)
		{
		  if(Sources->data[sn][Sources->datadict["primary_ID"]]==Sources->data[i][Sources->datadict["ID"]])
		    {
		      primarysn = i;
		      break;
		    }
		}

	      if(primarysn==-1)
		{
		  cout << "Problem finding primary star for multiple star where a companion was selected as the main source (sn,ID,primary_ID): ("
		       << sn << "," << Sources->data[sn][Sources->datadict["ID"]] << "," << Sources->data[sn][Sources->datadict["primary_ID"]] << ")" << endl;
		  exit(1);
		}

	      for(int i=primarysn;i<Sources->data.size();i++)
		{
		  if(i==primarysn || Sources->data[i][Sources->datadict["primary_ID"]]==Sources->data[sn][Sources->datadict["ID"]])
		    {
		      if(i!=sn)
			{
			  if(Paramfile->verbosity>1) cout << "Adding star " << i << " to source system." << endl;
			  Event->scompanions.push_back(i);
			  //add companion properties to the event data here?

			} //end if i!=sn
		    } //end if i==primarysn
		  else
		    {
		      break;
		    } //not sure I understand this break
		} //end for over sources
	    } //end else isbinary==1

	  //vector<double> scomp_rs;
	  //vector<double> scomp_s, scomp_alpha;
	  //vector<double> lcomp_s, lcomp_q;

	  
	} //end isbinary>0
      else
	{
	  //add dummy companion properties to the event data here?
	}
    }

  if(Paramfile->multiple_sources)
    {
      //add companion properties to the event data here?
      for(auto sc : Event->scompanions)
	{
	  Event->scomp_rs.push_back((Sources->data[sc][Sources->RADIUS] * Rsun / Sources->data[sc][Sources->DIST]) / Event->thE);
	  double P = pow(10,Sources->data[sc][Sources->datadict["combined_logP"]])/DAYINYR;
	  double M1 = Sources->data[sn][Sources->datadict["Mass"]];
	  double M2 = Sources->data[sc][Sources->datadict["Mass"]];
	  double acomb = pow(P*P*(M1+M2),1.0/3.0);
	  double a1 = M2/(M1+M2) * acomb;
	  double a2 = M1/(M1+M2) * acomb;
	  Event->scomp_s.push_back(acomb/(Event->thE * Sources->data[sn][Sources->DIST]));
	  Event->scomp_alpha.push_back(360.0*ran2(idum));
	  Event->scomp_phase.push_back(360.0*ran2(idum));
	  double rnd = ran2(idum);
	  Event->scomp_inc.push_back(180*(rnd<0.5?acos(2*rnd):-acos(2-2*rnd))/PI);
	  Event->scomp_fsofs1.push_back(vector<double>());
	  for(int filt=0; filt<=Paramfile->Nfilters;filt++)
	    {
	      double magnitude1 = Sources->mags[sn][filt];
	      double magnitude2 = Sources->mags[sc][filt];
	      double fs2ofs1 = pow(10,-0.4*(magnitude2-magnitude1));
	      Event->scomp_fsofs1.back().push_back(fs2ofs1);
	    }
	}
    }

  if(Paramfile->multiple_lenses)
    {
      
    }

  //The random parameters

  Event->t0 = double(Paramfile->NUM_SIM_DAYS)*ran2(idum);
  Event->t0range = double(Paramfile->NUM_SIM_DAYS);
  Event->weight_scale = 1.0;
  Event->alpha = 360.0 * ran2(idum);

  //u0 will be calculated after we know the blending

  //calculate the fundamental microlensing properties
    
  //double tE, rE, thE, piE, rs, mu, vt;

  //einstein radius
  //cout << x << endl;

  //in AU
  Event->rE = rEsun * sqrt(Lenses->data[ln][Lenses->MASS] 
			   * Sources->data[sn][Sources->DIST] * (1-x) * x);

  //in mas
  Event->thE = Event->rE/Lenses->data[ln][Lenses->DIST];

  //relative ls proper motion - lens motion relative to the source
  //in mas/yr
  //calculate the heliocentric relative proper motion
  pmgal[0] = Lenses->data[ln][Lenses->MUL]-Sources->data[sn][Sources->MUL];
  pmgal[1] = Lenses->data[ln][Lenses->MUB]-Sources->data[sn][Sources->MUB];

  //work out its absolute value
  Event->murel_l = pmgal[0];
  Event->murel_b = pmgal[1];
  Event->murel = qAdd(pmgal[0],pmgal[1]);

  //calculate the parallax
  Event->piE = (1-x)/Event->rE;
  Event->piEN = 0.0; //These will be computed if parallax is used
  Event->piEE = 0.0;



  

  //in km s-1
  Event->vt = Event->murel * Lenses->data[ln][Lenses->DIST] * AU/1000.0 / SECINYR;

  //in days
  Event->tE_h = DAYINYR * Event->thE / Event->murel;

  //in Einstein radii
  Event->rs = (Sources->data[sn][Sources->RADIUS] * Rsun / Sources->data[sn][Sources->DIST]) / Event->thE;
  //radius (Rsun) -> AU / Ds (kpc) -> mas / thetaE (mas) = ratio

  //rate weighting
  Event->raww = 2.0 * Event->thE * Event->murel;

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
  //Event->w = Event->u0max*Event->raww;
  
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

double inSeasont0range(struct filekeywords* Paramfile, struct obsfilekeywords World[],struct event* Event)
{
  //Compute the fractional coverage of dates where an event can happen
  double ngood=0;
  double nall=0;
  double t=0.5/4.0;
  long i=0;
  vector<long> maxi(Paramfile->numobservatories);
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++) maxi[obsidx]=World[obsidx].weatherSequence.size();
  while(t<Paramfile->NUM_SIM_DAYS)
	{
	  int thisgood=0;
	  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  if(World[obsidx].weatherSequence[i%maxi[obsidx]]>0) thisgood=1;
		}
	  ngood+=thisgood;
	  nall++;
	  t+=0.25;
	  i++;
	}
  //for(int i=0;i<World[obsidx].weatherSequence.size();i++)
  return double(ngood)/double(nall) * Paramfile->NUM_SIM_DAYS;
}


//void setupParallax(double tref, struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event *Event, struct slcat *Sources, struct slcat *Lenses)
void setupParallax(struct filekeywords* Paramfile, struct obsfilekeywords World[], struct event *Event, struct slcat *Sources, struct slcat *Lenses)
{
  if(Paramfile->verbosity>1) cout << __FUNCTION__ << endl;
  int sn = Event->source;
  int ln = Event->lens;
  double tref = Paramfile->tref;

  coords c;

  if(Paramfile->verbosity>0)
    {
      cout << "setupParallax" << endl;
    }

  Event->pllx.resize(Paramfile->numobservatories);
  
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //Prepare the epochs correctly
      //vector<double> jdepochs(1,Paramfile->simulation_zerotime);

      if(Paramfile->verbosity>3) Event->pllx[obsidx].debug=1;		  
      Event->pllx[obsidx].reset();
	  
      Event->pllx[obsidx].set_lb(Event->l, Event->b);
      Event->pllx[obsidx].setup_reference_frame(Paramfile->simulation_zerotime+tref,&World[0].orbit);
      Event->pllx[obsidx].set_orbit(&World[obsidx].orbit);
      //Event->pllx[obsidx].set_pm_lb(Lenses->data[ln][MUL]-Sources->data[sn][MUL],Lenses->data[ln][MUB]-Sources->data[sn][MUB]);
      //Event->pllx[obsidx].set_piE(Event->piE);
      //Event->pllx[obsidx].set_tE_h(Event->tE_h);
      Event->pllx[obsidx].provide_murel_h_lb(Lenses->data[ln][Lenses->MUL]-Sources->data[sn][Sources->MUL],
					     Lenses->data[ln][Lenses->MUB]-Sources->data[sn][Sources->MUB],
					     Event->piE, Event->thE);
	  
      //Do this for a dummy epoch now, but redo this at the end of 
      //timeSequencer
      //cout << "Observatory " << obsidx << " using " << World[obsidx].jd.size() << " epochs." << endl;
      if(Paramfile->verbosity>3) cout << "obsidx=" << obsidx << " jdtimes.size()=" << Event->jdtimes[obsidx].size() << endl;
      Event->pllx[obsidx].load_epochs(&Event->jdtimes[obsidx]);
      Event->pllx[obsidx].compute_NEshifts();
      Event->pllx[obsidx].compute_tushifts();
      //Event->pllx[obsidx].initialize();
	  
      if(Paramfile->verbosity>3)
	{
	  for(int idx=0;idx<int(Event->jdtimes[obsidx].size());idx++)
	    {
	      printf("NEShift event %d obs %d %14.6f %g %g\n",Event->id,obsidx,Event->jdtimes[obsidx][idx],Event->pllx[obsidx].Nshift[idx] + (Event->pllx[obsidx].epochs[idx]-Event->pllx[obsidx].tref)*Event->pllx[obsidx].vref[0],Event->pllx[obsidx].Eshift[idx] + (Event->pllx[obsidx].epochs[idx]-Event->pllx[obsidx].tref)*Event->pllx[obsidx].vref[1]);
	    }
	}

      
      if(Paramfile->verbosity>0)
	{
	  cerr << "Observatory " << obsidx << ":" << endl;
	  cerr << "Orbit: " << endl;
	  cerr << World[obsidx].orbit[0].xh << " " << World[obsidx].orbit[0].yh << " " << World[obsidx].orbit[0].zh << endl;
	  
	  cerr << "tref = " << Event->pllx[obsidx].tref << endl;
	  cerr << "l,b,a,d,ra,dec = " << Event->pllx[obsidx].l << " " << Event->pllx[obsidx].b << " " << Event->pllx[obsidx].a << " " << Event->pllx[obsidx].d << " " << Event->ra << " " << Event->dec << endl;
	  cerr << "helio_frame mul, mub, mua, mud, mulam, mubet = " << Event->pllx[obsidx].mul_h << " " << Event->pllx[obsidx].mub_h << " " << Event->pllx[obsidx].mua_h << " " << Event->pllx[obsidx].mud_h << " " << Event->pllx[obsidx].mulam_h << " " << Event->pllx[obsidx].mubet_h << " " << endl;
	  cerr << "refer_frame mul, mub, mua, mud, mulam, mubet = " << Event->pllx[obsidx].mul_r << " " << Event->pllx[obsidx].mub_r << " " << Event->pllx[obsidx].mua_r << " " << Event->pllx[obsidx].mud_r << " " << Event->pllx[obsidx].mulam_r << " " << Event->pllx[obsidx].mubet_r << " " << endl;
	  cerr << "piEN, piEE, piEll, piErp, piE = " << Event->pllx[obsidx].piEN << " " << Event->pllx[obsidx].piEE << " " << Event->pllx[obsidx].piEll << " " << Event->pllx[obsidx].piErp << " " << Event->pllx[obsidx].piE << endl;
	  cerr << "tE_h, tE_r = " << Event->pllx[obsidx].tE_h << " " << Event->pllx[obsidx].tE_r << endl;
	}
      
    }
  
  Event->piEN = Event->pllx[0].piEN; //*Event->piE;
  Event->piEE = Event->pllx[0].piEE; //*Event->piE;
  Event->tE_h = Event->pllx[0].tE_h;
  Event->tE_r = Event->pllx[0].tE_r;

  //Compute the event rate weighting
  Event->w = Event->raww * Event->weight_scale * Event->u0max * (Event->t0range/365.25) / (Event->tE_r/Event->tE_h);
  
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
  Event->obsgroupoutputheader.clear();
  Event->obsgroupoutputheader.resize(Event->obsgroups.size(),string(""));
}

