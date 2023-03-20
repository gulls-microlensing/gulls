#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "singleLens.h"

#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

#include<time.h>
#include<vector>

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
  double alpha, cosa, sina, xcom;
  double z1,z2;

  Event->Amax=-1;
  Event->umin=1e50;

  int idx,obsidx;

  int fsflag;
  int errflag;

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
  double tt, uu;

  fs.reset();

  Event->lcerror=0;
  errflag=0;

  vector<int> idxshift;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];

      if(Paramfile->identicalSequence && obsidx>0)
	{
	  //lightcurve is identical from observatory to observatory
	  amp = Event->Atrue[idx-idxshift[obsidx]];
	}
      else
	{

	  /* Compute the magnification */

	  //rotate the parallax shifts into the t,u plane
	  tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
	  uu = Event->u0;

	  if(Paramfile->pllxMultiplyer)
	    {
	      //tt += Event->pllx[obsidx].tshift(idx-idxshift[obsidx]);
	      //uu += Event->pllx[obsidx].ushift(idx-idxshift[obsidx]);
	      tt += Event->pllx[obsidx].tshift(Event->jdepoch[idx]);
	      uu += Event->pllx[obsidx].ushift(Event->jdepoch[idx]);
	    }

	  Event->umin=min(Event->umin,qAdd(tt,uu));

	  xsCenter = tt*cosa - uu*sina + xcom;
	  ysCenter = tt*sina + uu*cosa;
	  
	  u=qAdd(xsCenter-xcom,ysCenter);
	  if(fsflag) fsflag=0;

	  amp = fs.fsmag(cd(xsCenter,ysCenter), fsflag);
	  
	  if(fsflag)
	    {
	      fsflag=0;
	      magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, 
		       &eps, &fsflag); 
	    }
	  Event->lcerror = fsflag;
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
  
    }
}

