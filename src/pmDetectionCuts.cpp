#include "lightcurveFitter.h"
#include "detectionCuts.h"

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  int obsidx;

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  //Determine if an event is detected 

  //Fit a PS lightcurve
  Event->deterror=0;
  Event->detected=0;
  lightcurveFitter(Paramfile, Event);

  /* If finite source star effects important, fit FS lightcurve */
  if(abs(Event->u0) < 10*Event->rs && Event->PSPL.chisq>Paramfile->minChiSquared)
    {
      Event->flag_needFS=1;
      Event->deterror=0;
      Event->detected=0;
      lightcurveFitter_FS(Paramfile, Event); 
    }

  if(Event->deterror) 
    {
      fprintf(stderr,"\nDiscarding event %d (Failed fit FS, errval=%d) \n", 
	      Event->id, Event->deterror);
    }
  else
    {
      //did we detect it?
      if((!Event->flag_needFS && Event->PSPL.chisq>Paramfile->minChiSquared)
	 || (Event->flag_needFS && Event->FSPL.chisq>Paramfile->minChiSquared))
	Event->detected=1;
    }

  //Work out whether there are nearby stars with high fluxes 
  // - rebuild the images and look closely

  double bg,src,lens,stars,err;
  int nblends;
  vector<double> blendmag;
  vector<double> blendx, blendy;
  int Nsub, xsub, ysub;
  int satflag;
  double scalefac;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //reset the detectors
      World[obsidx].im.reset_image();
      World[obsidx].im.reset_detector();
      nblends=0;
      blendmag.clear();
      blendx.clear();
      blendy.clear();
      Nsub = World[obsidx].im.psf.Nsub;
      xsub = Event->xsub[obsidx];
      ysub = Event->ysub[obsidx];
      scalefac =  World[obsidx].im.psf.pixscale/Nsub/World[obsidx].im.fwhm;
      cout << obsidx << " " << World[obsidx].im.fwhm << " " << scalefac << endl;

      int filter = World[obsidx].filter;

      //Add the background
      World[obsidx].im.set_background(20.0-2.5*log10(World[obsidx].constbackground + World[obsidx].zodiflux[0]));
      World[obsidx].im.addbg();
      //photometry of just the background
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&bg, &err, &satflag, false);

      //Add the source
      World[obsidx].im.addstar(xsub, ysub,
			       Sources->mags[Event->source][filter]);
      //photometry of with the source
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&src, &err, &satflag, false);

      //Add the lens
      World[obsidx].im.addstar(xsub, ysub,
			       Lenses->mags[Event->lens][filter]);
      //photometry of with the lens
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&lens, &err, &satflag, false);

      //Add the background stars
      for(int i=0;i<Event->sl[obsidx].nstars;i++)
	{
	  World[obsidx].im.addstar(Event->sl[obsidx].x[i], 
				   Event->sl[obsidx].y[i], 
				   Event->sl[obsidx].mag[i]);
	  if(scalefac*qAdd(Event->sl[obsidx].x[i]-xsub,Event->sl[obsidx].y[i]-ysub)<=2)
	    {
	      nblends++;
	      blendmag.push_back(Event->sl[obsidx].mag[i]);
	      blendx.push_back(Event->sl[obsidx].x[i]/double(Nsub));
	      blendy.push_back(Event->sl[obsidx].y[i]/double(Nsub));
	    }
	}
      //photometry of with everything
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], 
					World[obsidx].mintexp, 1, 
					&stars, &err, &satflag, false);

      //remove the background
      World[obsidx].im.subbg();

      //remove the other contributions
      stars-=lens;
      lens-=src;
      src-=bg;

      //Work out the errors on photometry/astrometry of the host

      //store the data
      Event->data.push_back(xsub/double(Nsub));
      Event->data.push_back(ysub/double(Nsub));
      Event->data.push_back(bg);
      Event->data.push_back(src);
      Event->data.push_back(lens);
      Event->data.push_back(stars);
      Event->data.push_back(satflag);
      Event->data.push_back(nblends);
      for(int i=0;i<nblends;i++) 
	{
	  Event->data.push_back(blendx[i]);
	  Event->data.push_back(blendy[i]);
	  Event->data.push_back(blendmag[i]);
	}
    }
}
