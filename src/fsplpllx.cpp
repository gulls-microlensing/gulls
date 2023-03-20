#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "singleLens.h"
#include "lightcurveFitter.h"

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
  double slope;
  double tN;
  double theta, cosB, sinB, xcom, cosT, sinT;
  double z1,z2;

  Event->Amax=-1;

  int idx,obsidx;

  int fsflag;
  int errflag;

  double u0, du, dtN;
  double rot, crot, srot; //angle between E,N and ||,|_

  //work out the event parameters in the fortran parametrization

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];	    /* separation */
  z1 = -m1*a;
  z2 = (1-m1)*a;
  rs = Event->rs;	            /* source size */
  slope = Event->alpha*TO_RAD;	    /* slope of the trajectory */
  Gamma = Event->gamma;	            /* limb-darkening profile */

  cosB = cos(slope); sinB = sin(slope);
  theta = slope - PI/2.0;
  cosT = cos(theta); sinT = sin(theta);
  /*xcom = a*(1-2*m1);*/
  xcom = -m1*a;  /* Origin is primary lens*/

  double u;

  Event->lcerror=0;
  errflag=0;

  vector<int> idxshift;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  //work out the rotation for applying parallax
  rot = -atan2(Event->piEN,Event->piEE);
  crot = cos(rot); srot=sin(rot);

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  //Event->piE = 0;

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
	  du = (crot*Event->Eshift[idx] + srot*Event->Nshift[idx])*Event->piE;
	  dtN= (-srot*Event->Eshift[idx] + crot*Event->Nshift[idx])*Event->piE;

	  tN = (Event->epoch[idx] - Event->t0) / Event->tE + dtN;
	  u0 = Event->u0+du;

	  xsCenter = u0*cosB - tN*sinB + xcom;
	  ysCenter = u0*sinB + tN*cosB;
	  
	  u=qAdd(xsCenter-xcom,ysCenter);
	  if(fsflag) fsflag=0;

	  amp = wittFSMagnification(u,rs);
	}
	  
      Event->Atrue[idx] = amp;
 
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

      //keep track of highest magnification
      if(amp>Event->Amax) 
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}
  
    }
}

