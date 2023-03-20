#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "singleLens.h"
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

  vector<int> obsoffset(Paramfile->numobservatories,0);

  Event->Amax=-1;
  Event->umin=1e50;

  int idx,obsidx;

  int errflag;

  double tt, uu;

  //work out the event parameters in the fortran parametrization

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];	    /* separation */
  rs = Event->rs;	            /* source size */
  alpha = Event->alpha*TO_RAD;	    /* slope of the trajectory */
  Gamma = Event->gamma;	            /* limb-darkening profile */

  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  /*xcom = -m1*a;  /* Origin is primary lens*/
  xcom = -a; //Origin is the primary lens, calculation done in planet frame
  //xcom = a*(1-2*m1); /* Origin is the Center of mass */

  Event->lcerror=0;
  errflag=0;

  if(abs(Event->params[7])>1 || Event->params[9]==0)
    {
      Event->allsat=9;
    }

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat) 
    {
      return;
    }

  vector<int> idxshift;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

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

	  tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
	  uu = Event->u0;

	  if(Paramfile->pllxMultiplyer)
	    {
	      //if(tt<3)
		//cerr << tt << " " << uu << " " <<Event->pllx[obsidx].tshift(idx-idxshift[obsidx]) << " "  << Event->pllx[obsidx].ushift(idx-idxshift[obsidx]) << endl;
	      //tt += Event->pllx[obsidx].tshift(idx-idxshift[obsidx]);
	      //uu += Event->pllx[obsidx].ushift(idx-idxshift[obsidx]);
	      tt += Event->pllx[obsidx].tshift(Event->jdepoch[idx]);
	      uu += Event->pllx[obsidx].ushift(Event->jdepoch[idx]);
	    }

	  Event->umin=min(Event->umin,qAdd(tt,uu));

	  xsCenter = tt*cosa - uu*sina + xcom;
	  ysCenter = tt*sina + uu*cosa;

	  magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
		   &errflag); 
	  //xsCenter-=xcom; amp = pacAmp(qAdd(xsCenter,ysCenter));//for testing
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

  //if there has been an error - try the backup generator
  if(Event->lcerror)
    {
      Event->lcerror=0;
      backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
    }
}

