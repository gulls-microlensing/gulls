#include <ctime>
#include<time.h>

#include "structures.h"
#include "columnCodes.h"
#include "constdefs.h"
#include "strfns.h"
#include "info.h"

#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

void backupGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
  //A slower, but potentially more robust lightcurve generator for use when
  //the main generator fails

  char str[100];

  double m1, a;
  double xsCenter, ysCenter, rs, Gamma, z1, z2;
  double amp;
  double alpha, cosa, sina, xcom;

  double tt,uu;

  int idx,obsidx;
  
  int fsflag=0;
  int errflag=0;

  //set umin and Amax to unreasonable numbers?
  Event->umin=1e50;
  Event->Amax=-1;  

  //work out the event parameters in the fortran parametrization
  
  m1 = 1.0/(1 + Event->params[QQ]); /* mass of the first lens m1+m2=1 */
  a = Event->params[SS];	/* separation */
  // what are z1 and z2? Lens positions?
  z1 = -m1*a;
  z2 = (1-m1)*a;
  // gathered from event structure
  rs = Event->rs;	/* source size */
  alpha = Event->alpha*TO_RAD;	/* slope of the trajectory */
  Gamma = Event->gamma;	/* limb-darkening profile */

  cosa = cos(alpha); sina = sin(alpha);
  /*xcom = a*(1-2*m1);*/
  xcom = -m1*a;  /* Origin is primary lens*/

  //what is src_cld? in headers drectory, class for circ limb darkened source.
  src_cld src = src_cld(rs,0);
  //lens_binary is another class in header dir
  lens_binary lens = lens_binary(m1, z1, 1-m1, z2);
  //same thing, but finiteSource class is in fs.h
  finiteSource fs = finiteSource(&src,&lens,5.0e-4,1.0e-3);

  //must be some parameter that need to be reevaluated to see if fs is relevant
  fs.reset();

  Event->lcerror=0;
  errflag=0;

  //baseline may not be saturated, but all photometry may still be
  Event->allsat=1; 

  //set the stopwatch
  time_t starttime = time(NULL);
  time_t currtime=0;

  //Calculate the lightcurve
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];

      /* Compute the magnification */

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
      
      amp = fs.fsmag(cd(xsCenter,ysCenter), fsflag);

      if(fsflag) Event->lcerror = fsflag;

      Event->Atrue[idx] = amp;
 
      currtime = time(NULL);
      if(difftime(currtime,starttime)>MAX_TIME_PER_LC)
	{
	  Event->lcerror = 4999;
	  cerr << __FILE__ << " " << __FUNCTION__ << " Error: Lightcurve took too long to generate" << endl;
	  break;
	}

      if( Event->lcerror != 0)  
	{
	  Event->lcerror = 4001;
	  sprintf(str,"\nerror caught from fsmag errval:%d",Event->lcerror);
	  logfile_ptr << Event->lcerror << endl;
	  logfile_ptr << Event->u0 << " " << Event->tE_r << " " 
		      << Event->t0 << " " << Event->params[QQ] << " " 
		      << Event->params[SS] << " " << Event->rs << " " 
		      << xsCenter << " " << ysCenter << endl;
	  fmtline(str,WIDTH,"(lightcurveGenerator)");
	  errorHandler(Event->lcerror);
	  break;
	}

      //keep track of highest magnification
      if(amp>Event->Amax) 
	{
	  Event->Amax = amp;
	  Event->peakpoint = idx;
	}
  
    }

}
