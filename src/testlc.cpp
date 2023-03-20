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

int main()
{
  char str[100];

  double m1, a;
  double xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double alpha, cosa, sina, xcom;

  int idx,obsidx;

  int errflag;

  double tt, uu;

  double u0=-0.07;
  double t0 = 0.01;
  double tE = 1;

  //work out the event parameters in the fortran parametrization

  

  m1 = 1.0/(1 + 0.01); /* mass of the first lens m1+m2=1 */
  a = 1.05;	    /* separation */
  rs = 0.001;	            /* source size */
  alpha = 41*TO_RAD;	    /* slope of the trajectory */
  Gamma = 0;	            /* limb-darkening profile */

  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  xcom = -m1*a;  /* Origin is primary lens*/
  //xcom = a*(1-2*m1); /* Origin is the Center of mass */

  int lcerror=0;
  errflag=0;

  //Calculate the lightcurve
  for(double t=-1.5;t<1.5;t+=0.001)
    {
      /* Compute the magnification */

      tt = (t - t0) / tE;
      uu = u0;

      xsCenter = tt*cosa - uu*sina + xcom;
      ysCenter = tt*sina + uu*cosa;
      
      magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
	       &errflag); 
      //xsCenter-=xcom; amp = pacAmp(qAdd(xsCenter,ysCenter));//for testing
      cout << t << " " << amp << endl;
    }

  //if there has been an error - try the backup generator
  //  if(Event->lcerror)
  // {
  //  Event->lcerror=0;
  //  backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
  //}
}

