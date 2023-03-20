#include "lightcurveGenerator.h"
#include "astroFns.h"

#include<fstream>

#define DEBUGVAR 0

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, int *errval, ofstream& logfile_ptr, int lcid, char* eventprefix)
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

  double nc, err, baseline;
  double ampmag;
  int satflag;

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

  //cout << "input: " << Event->params[0] << " " << Event->params[1] << " " << Event->params[2] << " " << Event->params[3] << " " << Event->params[4] << " " << Event->params[5] << " " << Event->params[6] << endl;
  //cout << "mlparams: " << m1 << " " << a << " " << rs << " " << slope << " " << Gamma << " " << xcom << endl; 

  *errval=0;

  //if the event is saturated in each band, no need to calculate the lightcurve
  if(Event->nepochs==0||Event->allsat) 
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
      *errval = 0;

      obsidx = Event->obsidx[idx];

      /* Compute the magnification */
      tN = (Event->epoch[idx] - Event->t0) / Event->tE;

      xsCenter = Event->u0*cosB + tN*cosT + xcom;
      ysCenter = Event->u0*sinB + tN*sinT;

      magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
	       errval); 
 
      if( (*errval) != 0) 
	{
	  sprintf(str,"\nerror caught from magfunc_  errval:%d",*errval);
	  logfile_ptr << *errval << endl;
	  logfile_ptr << Event->u0 << " " << Event->tE << " " 
		      << Event->t0 << " " << Event->params[QQ] << " " 
		      << Event->params[SS] << " " << Event->rs << " " 
		      << xsCenter << " " << ysCenter << endl;
	  fmtline(str,WIDTH,"(lightcurveGenerator)");
	  errorHandler(*errval);
	  break;
	}

      filter = World[obsidx].filter;

      //add the background
      if(!CONSTANT_BACKGROUND) World[obsidx].im.addbg();

      ampmag = Sources->mags[sn][filter]-2.5*log10(amp);
      //add the star
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       ampmag);
      //photometry
      World[obsidx].im.ideal_photometry(Event->xpix[obsidx], 
					Event->ypix[obsidx], Event->texp[idx], 
					Event->nstack[idx], &nc, &err, 
					&satflag, bool(CONSTANT_BACKGROUND));
      //cout << obsidx << " " << amp << " " << ampmag << " " << Sources->mags[sn][filter] << " " << nc / (Event->texp[idx]*Event->nstack[idx]) << endl;
      if(Event->allsat && !satflag) Event->allsat = 0;
      if(Event->allsatobs[obsidx] && !satflag) Event->allsatobs[obsidx] = 0;
      //subtract the star
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       ampmag);
      //subtract the background
      if(!CONSTANT_BACKGROUND) World[obsidx].im.subbg();

      //store the results
      baseline = Event->baselineFlux[obsidx] * Event->texp[idx] 
	* Event->nstack[idx];
      Event->Atrue[idx] = nc/baseline;
      Event->Atrueerr[idx] = err/baseline;

      if(Paramfile->outputLightcurve)
	{
	  //reperform the photometry, but with noise
	  //add the background
	  if(!CONSTANT_BACKGROUND) World[obsidx].im.addbg();

	  //add the star
	  World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);
	  //photometry
	  World[obsidx].im.quick_photometry(Event->xpix[obsidx], 
					    Event->ypix[obsidx], 
					    Event->texp[idx], 
					    Event->nstack[idx], &nc, &err, 
					    &satflag);
	  Event->Aobs[idx] = nc/baseline;
	  Event->Aerr[idx] = err/baseline;
	  Event->nosat[idx] = !satflag;

	  World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   ampmag);
	  //subtract the background
	  if(!CONSTANT_BACKGROUND) World[obsidx].im.subbg();
	}
      else
	{
	  Event->Aobs[idx] = Event->Atrue[idx];
	  Event->Aerr[idx] = Event->Atrueerr[idx];
	  Event->nosat[idx] = !satflag; //nosat is the oposite of satflag
	}

      //keep track of highest magnification
      if(amp>Event->Amax) Event->Amax = amp;
  
    }
}

