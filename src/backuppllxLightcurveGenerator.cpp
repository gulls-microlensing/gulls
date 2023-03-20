#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include<time.h>

#include<fstream>

#define DEBUGVAR 0

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[100];

  double m1, a;
  double xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double slope;
  double tN;
  double theta, cosB, sinB, xcom, cosT, sinT;

  Event->Amax=-1;

  int idx,obsidx;
  
  int filter;
  int sn = Event->source;

  double baseline;
  double ampmag;
  int satflag;
  int errflag;

  timespec lcstart, lcend, photstart, photend;
  long nsec;
  double nci, ncs, erri, errs;  

  double u0, du, dtN;
  double rot, crot, srot; //angle between E,N and ||,|_

  //work out the event parameters in the fortran parametrization

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];	/* separation */
  rs = Event->rs;	/* source size */
  slope = Event->alpha*TO_RAD;	/* slope of the trajectory */
  Gamma = Event->gamma;	/* limb-darkening profile */

  cosB = cos(slope); sinB = sin(slope);
  theta = slope - PI/2.0;
  cosT = cos(theta); sinT = sin(theta);
  /*xcom = a*(1-2*m1);*/
  xcom = -m1*a;  /* Origin is primary lens*/

  Event->lcerror=0;
  errflag=0;

  //work out the rotation for applying parallax
  rot = -atan2(Event->piEN,Event->piEE);
  crot = cos(rot); srot=sin(rot);

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  //baseline may be unsaturated, but all photometry may still be
  Event->allsat=1;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++) 
    Event->allsatobs[obsidx]=1;

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];

      /* Compute the magnification */

      //rotate the parallax shifts into the t,u plane
      du = (crot*Event->Eshift[idx] + srot*Event->Nshift[idx])*Event->piE;
      dtN = (-srot*Event->Eshift[idx] + crot*Event->Nshift[idx])*Event->piE;

      tN = (Event->epoch[idx] - Event->t0) / Event->tE + dtN;
      u0 = Event->u0+du;

      xsCenter = u0*cosB + tN*cosT + xcom;
      ysCenter = u0*sinB + tN*sinT;

      //time the magnification calculation
      clock_gettime(CLOCK_REALTIME,&lcstart);

      magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
	       &errflag); 

      clock_gettime(CLOCK_REALTIME,&lcend);
      nsec = lcend.tv_nsec - lcstart.tv_nsec;
      Paramfile->lctime += double((lcend.tv_sec - lcstart.tv_sec) 
				  - (nsec<0?1:0))
	+ double(nsec<0?nsec+1000000000:nsec)*1.0e-9;
 
      if( errflag != 0) 
	{
	  sprintf(str,"\nerror caught from magfunc_  errval:%d", 
		  Event->lcerror);
	  logfile_ptr << Event->lcerror << endl;
	  logfile_ptr << Event->u0 << " " << Event->tE << " " 
		      << Event->t0 << " " << Event->params[QQ] << " " 
		      << Event->params[SS] << " " << Event->rs << " " 
		      << xsCenter << " " << ysCenter << endl;
	  fmtline(str,WIDTH,"(lightcurveGenerator)");
	  errorHandler(errflag);
	  Event->lcerror=errflag;
	  //break;
	}

      filter = World[obsidx].filter;

      //time the photometry
      clock_gettime(CLOCK_REALTIME,&photstart);

      //add the background
      //if(!CONSTANT_BACKGROUND) World[obsidx].im.addbg();
      World[obsidx].im.set_background(Event->backmag[idx]);
      World[obsidx].im.addbg();

      ampmag = Sources->mags[sn][filter]-2.5*log10(amp);
      //add the star
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       ampmag);
      //photometry
      World[obsidx].im.ideal_and_scattered_photometry(Event->xpix[obsidx],
						      Event->ypix[obsidx], 
						      Event->texp[idx],
						      Event->nstack[idx], 
						      &nci, &erri, &ncs, &errs,
						      &satflag);

      if(Event->allsat && !satflag) Event->allsat = 0;
      if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;

      //store the results
      baseline = Event->baselineFlux[obsidx] * Event->texp[idx] 
	* Event->nstack[idx];
      Event->Atrue[idx] = nci/baseline;
      Event->Atrueerr[idx] = erri/baseline;

      if(!Paramfile->idealPhotometry)
	{
	  //Use the noisy photometry
	  Event->Aobs[idx] = ncs/baseline;
	  Event->Aerr[idx] = errs/baseline;
	  Event->nosat[idx] = !satflag;
	}
      else
	{
	  Event->Aobs[idx] = Event->Atrue[idx];
	  Event->Aerr[idx] = Event->Atrueerr[idx];
	  Event->nosat[idx] = !satflag; //nosat is the oposite of satflag
	}

      //subtract the star
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       ampmag);
      //subtract the background
      World[obsidx].im.subbg();

      clock_gettime(CLOCK_REALTIME,&photend);
      nsec = photend.tv_nsec - photstart.tv_nsec;
      Paramfile->phottime += double((photend.tv_sec - photstart.tv_sec) 
				  - (nsec<0?1:0))
	+ double(nsec<0?nsec+1000000000:nsec)*1.0e-9;

      //keep track of highest magnification
      if(amp>Event->Amax) 
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}
  
    }

  //if there has been an error - try the backup generator
  if(Event->lcerror)
    {
      Event->lcerror=0;
      backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
    }
}

