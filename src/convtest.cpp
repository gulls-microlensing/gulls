#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"

#include<fstream>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

#define DEBUGVAR 1

extern "C"
{
  void magfunc_(double *m1, double *a, double *xsCenter,  double *ysCenter, double *rs, double *Gamma, double *amp, double *eps, int *errflag);
}

void fisherMatrix(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{
  char str[512];

  double m1, a;
  double xsCenter, ysCenter, rs, Gamma;
  double amp, eps=1.0e-3;
  double slope;
  double tN;
  double theta, cosB, sinB, xcom, cosT, sinT;
  double t0, tE, u0;
  double logq, logs, logtE;

  int idx,obsidx;
  
  int filter;
  int sn = Event->source;

  double nc, err, baseline;
  double ampmag;
  int satflag;
  int errflag;

  FILE* fmfile_ptr;
  FILE* ctfile_ptr;

  if(DEBUGVAR) cout << "Calculating fisher matrix for event number " << Event->id << endl;

  int Nparams=7;
  int Nobsparams = 2*Paramfile->numobservatories;
  int Ntotparams = Nparams + Nobsparams;
  // 0 = t0  1 = tE  2 = u0  3 = alpha  4 = s  5 = q  6 = rs  7 = F0  8 = fs
  double* step = new double[Ntotparams];
  double* F0 = new double[Paramfile->numobservatories];
  double* fs = new double[Paramfile->numobservatories];
  double A[2];
  double df;
  int dshift;
  double cstep;

  step[0] = 0.023; step[1] = 0.0023; step[2] = 0.0023; step[3] = 0.023;
  step[4] = 0.13; step[5] = 0.0023; step[6] = 0.0023;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      step[Nparams+2*obsidx] = 0.0013;
      step[Nparams+2*obsidx+1] = 0.0013;
    }

  for(int param=0;param<Ntotparams;param++)
    {
      if(DEBUGVAR) cout << "Parameter: " << param << endl;

      Event->lcerror=0;
      errflag=0;

      char ctfname[1000];

      sprintf(ctfname, "%s%s_%d_%d.det.fm.%d", Paramfile->outputdir,
	      Paramfile->run_name, Event->instance, Event->id, param);

  
      ctfile_ptr = fopen(ctfname,"w");

      //Calculate the lightcurve
      for(idx=0;idx<Event->nepochs;idx++)
	{

	  //work out the event parameters in the fortran parametrization
      
	  //mass of the first lens m1+m2=1
	  
	  m1 = 1.0/(1 + Event->params[QQ]); 
	  a = Event->params[SS];	//separation
	  rs = Event->rs;	        //source size
	  /* slope of the trajectory */
	  slope = Event->alpha*TO_RAD;
	  Gamma = Event->gamma;	/* limb-darkening profile */

	  cosB = cos(slope); sinB = sin(slope);
	  theta = slope - PI/2.0;
	  cosT = cos(theta); sinT = sin(theta);
	  /*xcom = a*(1-2*m1);*/
	  xcom = -m1*a;  /* Origin is primary lens*/
	  
	  t0 = Event->t0;
	  tE = Event->tE;
	  u0 = Event->u0;

	  for(int o=0;o<Paramfile->numobservatories;o++)
	    {
	      F0[o] = 1;
	      fs[o] = Event->fs[o];
	    }
	  
	  obsidx = Event->obsidx[idx];

	  //reset the parameters at each epoch - TO DO!!!!!!!!!!!!!!!!!!!!
	  
	  /* Compute the magnification */
	  tN = (Event->epoch[idx] - t0) / tE;
	  
	  xsCenter = u0*cosB + tN*cosT + xcom;
	  ysCenter = u0*sinB + tN*sinT;
	  
	  magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, &amp, &eps,
		   &errflag);

	  double u = sqrt(u0*u0+tN*tN);

	  if(idx%100==0 || (abs(amp-(u*u+2)/(u*sqrt(u*u+4)))/Event->Aerr[idx]>0.1&&!errflag))
	    {
	      for(int s=-144;s<=56;s++)
		{
		  cstep = step[param]*pow(2.0,double(s)/4.0);
		  if(param==4&&(cstep<1e-7||cstep>0.1)) continue;
		  else if(param>7&&param%2==0&&cstep>0.1) continue;
		  fprintf(ctfile_ptr,"%6g %6g ",Event->epoch[idx],cstep);


		  for(int delta=-1;delta<=1;delta+=2)
		    {


		      //mass of the first lens m1+m2=1
		      logq = log10(Event->params[QQ]) + (param==5?delta*cstep:0);
		      logs = log10(Event->params[SS]) + (param==4?delta*cstep:0);
		      logtE = log10(Event->tE) + (param==1?delta*cstep:0);
		      
		      m1 = 1.0/(1 + pow(10,logq)); 
		      a = pow(10,logs);	//separation
		      rs = pow(10,log10(Event->rs) + (param==6?delta*cstep:0));	        //source size
		      /* slope of the trajectory */
		      slope = Event->alpha*TO_RAD + (param==3?delta*cstep:0); 
		      Gamma = Event->gamma;	/* limb-darkening profile */

		      cosB = cos(slope); sinB = sin(slope);
		      theta = slope - PI/2.0;
		      cosT = cos(theta); sinT = sin(theta);
		      /*xcom = a*(1-2*m1);*/
		      xcom = -m1*a;  /* Origin is primary lens*/

		      t0 = Event->t0 + (param==0?delta*cstep:0);
		      tE = pow(10,logtE);
		      if(abs(Event->u0)<1.0e-5) u0 = Event->u0 + (param==2?delta*cstep:0);
		      else u0 = pow(10,log10(Event->u0) + (param==2?delta*cstep:0));

		      for(int o=0;o<Paramfile->numobservatories;o++)
			{
			  F0[o] = 1 + (param==Nparams+2*o?delta*cstep:0);
			  fs[o] = pow(10, log10(Event->fs[o]) + (param==Nparams+2*o+1?delta*cstep:0));
			}

		      /* Compute the magnification */
		      tN = (Event->epoch[idx] - t0) / tE;
	  
		      xsCenter = u0*cosB + tN*cosT + xcom;
		      ysCenter = u0*sinB + tN*sinT;

		      //if(s==-144||s==-143)
		      //	printf("%d %f %f %f %f %f %f\n",s,m1,a,xsCenter,ysCenter,rs,Gamma);
	  
		      magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, 
			       &amp, &eps, &errflag);
		      if(s==-144) magfunc_(&m1, &a, &xsCenter, &ysCenter, &rs, &Gamma, 
			       &amp, &eps, &errflag);
		      A[(delta+1)/2] = F0[obsidx] + fs[obsidx]*(amp-1);
		    } //end for delta

		  //compute the numerical derivatives
		  df = 0.5*(A[1]-A[0])/cstep;

		  if(!errflag)
		    fprintf(ctfile_ptr,"%6g\n",df);
		  else
		    fprintf(ctfile_ptr,"(1/0)\n",df);
		  errflag=0;

		} //end for s
	      fprintf(ctfile_ptr,"\n");
	    } //end if idx%100==0 
 
	  if( errflag != 0) 
	    {
	      //sprintf(str,"\nerror caught from magfunc_  errval:%d", 
	      //	      Event->lcerror);
	      //fmtline(str,WIDTH,"(fisher)");
	      //errorHandler(errflag);
	      Event->lcerror=errflag;
	      errflag=0;
	      Event->lcerror=0;
	      //break;
	    }
	  
	} //end for idx
    }//nd for param
  
  delete[] step;
  delete[] F0;
  delete[] fs;

  //if there has been an error - try the backup generator
  //  if(Event->lcerror)
  //  {
  //    Event->lcerror=0;
  //    backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
  //  }
}

