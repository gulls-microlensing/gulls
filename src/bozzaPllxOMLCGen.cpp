#include "lightcurveGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "constdefs.h"
#include<time.h>
#include<vector>
#include<iomanip>
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
  double alpha, cosa, sina;
  
  vector<int> obsoffset(Paramfile->numobservatories,0);

  Event->Amax=-1;
  Event->umin=1e50;
  double lim_gamma=Event->gamma;
  double lcgen=Paramfile->LC_GEN;
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
  double zerotime = Event->t0;
  if(Paramfile->parameterization==1) zerotime = Event->tcroin;

  double x20 = a2 * cos(phase0);
  double y20 = a2 * sin(phase0) * cosinc;
  double s0 = Event->params[SS];
  double xCoM = s0*q/(1+q); //Position of the center of mass in the frame of the primary lens at t0 or tcroin

  
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
  Event->vbm->a1 = lim_gamma;             /*  Linear limb-darkening coefficient.*/
 


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
  if(Paramfile->verbosity>=4)
    {
            Event->vbm_rootaccuracy.resize(Event->nepochs);
            Event->vbm_squarecheck.resize(Event->nepochs);
            Event->vbm_therr.resize(Event->nepochs);
    }  
  vector<int> idxshift;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  int shiftedidx;

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      shiftedidx = idx-idxshift[obsidx];

      if(Paramfile->identicalSequence && obsidx>0)
	{
	  //lightcurve is identical from observatory to observatory
	  amp = Event->Atrue[shiftedidx];
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
	      tt += Event->pllx[obsidx].tshift[shiftedidx];
	      uu += Event->pllx[obsidx].ushift[shiftedidx];
	    }

	  Event->umin=min(Event->umin,qAdd(tt,uu));
	  

	  //orbital calculations
	  phase = phase0 + 2*pi*(Event->epoch[idx] - zerotime) / period;
	  x1 = a1 * cos(phase+pi);
	  y1 = a1 * sin(phase+pi) * cosinc;
	  x2 = a2 * cos(phase);
	  y2 = a2 * sin(phase) * cosinc;
	  Event->xl1[idx] = x1;	  Event->yl1[idx] = y1;
	  Event->xl2[idx] = x2;	  Event->yl2[idx] = y2;
	  s = qAdd(x2-x1,y2-y1);
	  rot = atan2(y2,x2) - atan2(y20,x20);
	  cosrot = cos(-rot); sinrot = sin(-rot);

	  //source position
	  xsin = tt*cosa - uu*sina - xCoM; //as viewed from earth in non-rotating frame
	  ysin = tt*sina + uu*cosa;
	  Event->xs[idx] = xsin; Event->ys[idx] = ysin;
	  xsrot = xsin*cosrot - ysin*sinrot; //in frame rotating with binary
	  ysrot = xsin*sinrot + ysin*cosrot;

	  //VBB uses CoM as the origin
	  //amp = VBBL.BinaryMagDark(s,q,xsrot,ysrot,rs,Gamma,eps);
          //amp = VBBL.BinaryMag2(s, q, xsrot,ysrot,rs);
	  if(Paramfile->verbosity>=3)
	    cout << setprecision(16) << s << " " << q << " " << xsrot
		 << " " << ysrot << " " << rs << setprecision(6)
		 << endl;
	  
	  amp = Event->vbm->BinaryMag2(s,q,xsrot,ysrot,rs);
          if(Paramfile->verbosity>=4)
            {
                 Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
                 Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
                 Event->vbm_therr[idx] = Event->vbm->therr;
             }
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
		      << xsrot << " " << ysrot << endl;
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

