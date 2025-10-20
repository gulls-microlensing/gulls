#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "singleLens.h"
#include<time.h>
#include<vector>

#include<fstream>
#include <sstream>
#include<iomanip>
#include <sys/stat.h>
#define DEBUGVAR 0

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  char str[512];

  double m1, a;
  double q;
  double xsCoM,  xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina,VBM_origin,Mao_origin ;
  int useVBB=1;
  vector<int> obsoffset(Paramfile->numobservatories,0);
  Event->Amax=-1;
  Event->umin=1e50;
  double ampoldlc, ampvbm, dif_over_amp;
  int idx,obsidx;
  double lim_gamma=Paramfile->LD_GAMMA;
  int errflag;

  double tt, uu;
  double lcgen=Paramfile->LC_GEN;
  if(DEBUGVAR) cout << "LC_GEN: " << lcgen << endl;

  Event->vbm->SetMethod(VBMicrolensing::Method::Multipoly);

  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];            /* separation */
  rs = Event->rs;                   /* source size */
  alpha = Event->alpha*TO_RAD;      /* slope of the trajectory */
  Gamma = Event->gamma;             /* limb-darkening profile */
  Event->vbm->a1 = lim_gamma;              /*  Linear limb-darkening coefficient.*/
  q = Event->params[QQ];
  cosa = cos(alpha); sina = sin(alpha);


  Mao_origin = -m1*a; //Translate from L1 origin to m1z2+m2z1=0 origin
  VBM_origin = (1-m1)*(-a); //Translate from L1 origin to center of mass origin

  fstream fstr;
  fstr.open("test.txt",fstream::out);
  
  
  int nn=4;
  
  double pr[] = {     //parameters
    0.0, 0.0, 1.0,    // First lens: x1_1, x1_2, m1
    1.0, -0.7, 1e-4,  // Second lens: x2_2, x2_2, m2
    2.0, 0.7, 1e-4,   // Third lens: x3_re, x3_im, m3
    0.6, -0.6, 1e-6   // Fourth lens: x4_re, x4_im, m4
  };


  if(Paramfile->verbosity>2) cout << "Set geometry" << endl;
  Event->vbm->SetLensGeometry(nn,pr);
  


  Event->lcerror=0;
  errflag=0;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0 || Event->allsat)
    {
      return;
    }

  vector<int> idxshift;
  int shiftedidx;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      shiftedidx=idx-idxshift[obsidx];

      if(Paramfile->identicalSequence && obsidx>0)

        {
          //lightcurve is identical from observatory to observatory
          amp = Event->Atrue[shiftedidx];
        }
      else
        {

          /* Compute the magnification */

          tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
          uu = Event->u0;

          if(Paramfile->pllxMultiplyer)
            {
              //tt += Event->pllx[obsidx].tshift(Event->jdepoch[idx]);
              //uu += Event->pllx[obsidx].ushift(Event->jdepoch[idx]);
              tt += Event->pllx[obsidx].tshift[shiftedidx];
              uu += Event->pllx[obsidx].ushift[shiftedidx];
            }


          Event->umin=min(Event->umin,qAdd(tt,uu));
          
          xsCoM = tt*cosa - uu*sina + VBM_origin; //coordinate shift to center of mass
          xsCenter = tt*cosa - uu*sina + Mao_origin;// coordinate shif to primary lens
          ysCenter = tt*sina + uu*cosa;
	  if(Paramfile->verbosity>3) cout << hexfloat << Event->epoch[idx] << " " << Event->t0 << " " << Event->tE_r << " " << tt << " " << xsCoM << " " << ysCenter << " " << rs << endl;
	  if(Paramfile->verbosity>3) fstr << hexfloat << Event->epoch[idx] << " " << Event->t0 << " " << Event->tE_r << " " << tt << " " << xsCoM << " " << ysCenter << " " << rs << endl;
	  amp = Event->vbm->MultiMag2(xsCoM, ysCenter,rs);
	  if(Paramfile->verbosity>3) cout << "Event->vbm->MultiMag2(xsCoM, ysCenter,rs); done" << endl;
	  Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
	  Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
	  Event->vbm_therr[idx] = Event->vbm->therr;
	  
	  if(Paramfile->verbosity>3) cout << amp << endl;
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

  fstr.close();

}
