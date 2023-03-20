#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "singleLens.h"

#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

#include<time.h>

#include<fstream>

#define DEBUGVAR 0

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[100];

  double m1, a;
  double xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina, xcom;
  double z1,z2;

  Event->Amax=-1;

  int idx,obsidx;

  int fsflag;
  int errflag;

  double uu, tt;

  //work out the event parameters in the fortran parametrization

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];	    /* separation */
  z1 = -m1*a;
  z2 = (1-m1)*a;
  rs = Event->rs;	            /* source size */
  alpha = Event->alpha*TO_RAD;	    /* slope of the trajectory */
  Gamma = Event->gamma;	            /* limb-darkening profile */

  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  xcom = -m1*a;  /* Origin is primary lens*/

  src_cld src = src_cld(rs,0);
  lens_binary lens = lens_binary(m1, z1, 1-m1, z2);
  finiteSource fs = finiteSource(&src,&lens,5.0e-4,eps);

  double u;

  fs.reset();

  Event->lcerror=0;
  errflag=0;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  //Calculate the lightcurve

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      double ldu = double(obsidx)*0.3;
      src.setgamma(2.0/3.0 / (1-1/ldu)); //no reset needed - source is same size
      for(idx=Event->nepochsvec[obsidx];idx<Event->nepochsvec[obsidx+1];idx++)
	{
	  obsidx = Event->obsidx[idx];

	  /* Compute the magnification */

	  tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
	  uu = Event->u0;

	  xsCenter = tt*cosa - uu*sina + xcom;
	  ysCenter = tt*sina + uu*cosa;

	  u=qAdd(xsCenter-xcom,ysCenter);
	  if(u>10)
	    {
	      //if u is large - assume we have a single lens
	      amp = pacAmp(u);
	    }
	  else
	    {
	      if(fsflag) fsflag=0;
	      
	      amp = fs.fsmag(cd(xsCenter,ysCenter), fsflag);
	      if(fsflag) Event->lcerror = fsflag;
	    }

	  Event->Atrue[idx] = amp;
 
	  if( errflag != 0) 
	    {
	      sprintf(str,"\nerror caught from magfunc_  errval:%d", 
		      Event->lcerror);
	      logfile_ptr << Event->lcerror << endl;
	      logfile_ptr << Event->u0 << " " << Event->tE_r << " " 
			  << Event->t0 << " " << Event->params[QQ] << " " 
			  << Event->params[SS] << " " << Event->rs << " " 
			  << xsCenter << " " << ysCenter << endl;
	      fmtline(str,WIDTH,"(lightcurveGenerator)");
	      errorHandler(errflag);
	      Event->lcerror=errflag;
	      //break;
	    }

	  //keep track of highest magnification
	  if(amp>Event->Amax) 
	    {
	      Event->Amax = amp;
	      Event->peakpoint = idx;
	    }
	} //end for idx
    } //end for observatory
}

