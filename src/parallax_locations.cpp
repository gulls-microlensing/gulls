#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "VBBinaryLensingLibrary.h"
#include "constdefs.h"
#include<time.h>
#include<vector>

#include<fstream>

#define DEBUGVAR 0

//extern "C"
//{
//  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
//}

//#include "singleLens.h"

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[100];

  double xsCenter, ysCenter, rs, Gamma=0.4;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina, xcom;
  int useVBB=1;
  VBBinaryLensing VBBL;

  vector<int> obsoffset(Paramfile->numobservatories,0);

  Event->Amax=-1;
  Event->umin=1e50;

  int idx,obsidx;

  int errflag;

  double tt, uu;

  double aorE = Event->params[AA]/Event->rE;
  double cosinc = cos(Event->params[INC]*TO_RAD);
  double phase0 = Event->params[PHASE]*TO_RAD;
  double period = Event->params[TT]*DAYINYR;
  double q=Event->params[QQ];
  double s;
  double a1 = q/(1+q)*aorE; //in rE
  double a2 = aorE - a1;
  
  double phase;
  double x1, y1; //host position in inertial frame in rE
  double x2, y2; //planet position in inertial frame in rE
  double rot, cosrot, sinrot; //rotation angle to subtract to put source in rotating frame
  double xsin, ysin; //source position in the inertial frame
  double xsrot, ysrot; //source position in  the rotating frame
  

  //work out the event parameters in the fortran parametrization

  rs = Event->rs;	            /* source size */
  alpha = Event->alpha*TO_RAD;	    /* slope of the trajectory */
  Gamma = Event->gamma;	            /* limb-darkening profile */

  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  //xcom = -m1*a;  /* Origin is primary lens*/
  //xcom = a*(1-2*m1); /* Origin is the Center of mass */

  Event->lcerror=0;
  errflag=0;

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

	  //parallax shifts in the fixed reference frame
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
	  

	  //orbital calculations
	  phase = phase0 + 2*pi*(Event->epoch[idx] - Event->t0) / period;
	  x1 = a1 * cos(phase+pi);
	  y1 = a1 * sin(phase+pi) * cosinc;
	  x2 = a2 * cos(phase);
	  y2 = a2 * sin(phase) * cosinc;

	  //Spacecraft orbit
	  double xeq=0.0; double xecl=0.0;
	  double yeq=0.0; double yecl=0.0;
	  double zeq=0.0; double zecl=0.0;
	  for(int jj=0;jj<(*Event->pllx[obsidx].orbit).size();jj++)
	    {
	      (*Event->pllx[obsidx].orbit)[jj].settime(Event->pllx[obsidx].epochs[Event->jdepoch[idx]]);
	      (*Event->pllx[obsidx].orbit)[jj].compute();
	      xeq += (*Event->pllx[obsidx].orbit)[jj].xeq;
	      yeq += (*Event->pllx[obsidx].orbit)[jj].yeq;
	      zeq += (*Event->pllx[obsidx].orbit)[jj].zeq;
	      xecl += (*Event->pllx[obsidx].orbit)[jj].xecl;
	      yecl += (*Event->pllx[obsidx].orbit)[jj].yecl;
	      zecl += (*Event->pllx[obsidx].orbit)[jj].zecl;
	    }
	  
	  Event->xl1[idx] = xeq;
	  Event->yl1[idx] = yeq;
	  Event->xl2[idx] = zeq;
	  Event->yl2[idx] = xecl;
	  Event->xs[idx] = yecl;
	  Event->ys[idx] = zecl;
	  s = qAdd(x2-x1,y2-y1);
	  rot = atan2(y2,x2);
	  cosrot = cos(-rot); sinrot = sin(-rot);

	  //source position
	  xsin = tt*cosa - uu*sina; //as viewed from earth in non-rotating frame
	  ysin = tt*sina + uu*cosa;
	  //Event->xs[idx] = xsin; Event->ys[idx] = ysin;
	  xsrot = xsin*cosrot - ysin*sinrot; //in frame rotating with binary
	  ysrot = xsin*sinrot + ysin*cosrot;

	  //VBB uses CoM as the origin
	  amp = VBBL.BinaryMagDark(s,q,xsrot,ysrot,rs,Gamma,eps);
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

}

