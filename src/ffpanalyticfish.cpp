#include "lightcurveGenerator.h"
#include "lightcurveFitter.h"
#include "astroFns.h"
#include "integerPowers.h"
#include "wittFSPL.h"

#include<fstream>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

#define DEBUGVAR 1

void fisherMatrix(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{
  double rs, Gamma;
  double amp;
  double t0, tE, u0;
  double tt,uu;
  double logtE=0, logrs=0;
  double pllx=Paramfile->pllxMultiplyer;
  vector<parallax> parlx;

  int idx,obsidx;
  
  double piEN=0, piEE=0, piE;

  FILE* fmfile_ptr;
  ofstream lc("lc.txt");
  lc.precision(16);

  if(DEBUGVAR) cout << "Calculating fisher matrix for event number " << Event->id << endl;

  int Nparams=5;
  if(pllx==0) Nparams=3;
  int Nobsparams = 2*Paramfile->numobservatories;
  int Ntotparams = Nparams + Nobsparams;
  // 0 = t0  1 = tE  2 = u0  3 = rs  4 = piEN 5 = piEE 6 = F0  7 = fs
  double* step = new double[Ntotparams];
  double* dF = new double[Ntotparams*Event->nepochs];
  double* fs = new double[Paramfile->numobservatories];
  int dshift;
  double dmult;
  double u;

  int pshift;

  vector<int> idxshift;
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    idxshift.push_back(Event->nepochsvec[obsidx]);

  step[0] = 1e-6; step[1] = 1e-6; step[2] = 1e-6; //step[3] = 1e-6;
  if(pllx!=0)
    {
      step[3] = 1e-6;
      step[4] = 1e-6;
    }
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      step[Nparams+2*obsidx] = 1e-6;
      step[Nparams+2*obsidx+1] = 1e-6;
      parlx.push_back(Event->pllx[obsidx]);
    }

  //first initialize the derivatives to zero so that we can accumulate the sum
  for(int param=0; param<Ntotparams; param++)
    {
      int pshift = param*Event->nepochs;
      for(idx=0; idx<Event->nepochs; idx++) dF[pshift+idx]=0;
    }

  cerr.precision(16);


  //Calculate derivatives numerically for the non-linear parameters
  for(int param=0;param<Nparams;param++)
    {
      pshift = param*Event->nepochs;

      if(DEBUGVAR) cout << "Parameter: " << param << endl;
      for(int delta=-1;delta<=+1;delta+=2)
	{
	  
	  dshift = Event->nepochs*(delta+1)/2;
	  dmult = 0.5*delta/step[param];

	  t0 = Event->t0 + (param==0?delta*step[0]:0);
	  logtE = log10(Event->tE_r) + (param==1?delta*step[1]:0);
	  tE = pow(10,logtE);
	  u0 = Event->u0 + (param==2?delta*step[2]:0);
	  //logrs = log10(Event->rs) + (param==3?delta*step[3]:0);
	  rs = Event->rs; //pow(10,logrs);
	  Gamma = Event->gamma;	/* limb-darkening profile */
	  if(pllx)
	    {
	      piEN = Event->piEN + (param==4?delta*step[4]:0);
	      piEE = Event->piEE + (param==5?delta*step[5]:0);
	      piE = qAdd(piEN,piEE);
	    }

	  cerr << param << " " << delta << " " << (delta>0?" ":"") << t0 << " " << tE << " " << u0 << " " << rs << " " << piEN << " " << piEE << endl;


	  //only adjust parallax parameters if necessary
	  if(pllx && (param==3 || param==4 || param == 1))
	    {
	      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  parlx[obsidx].set_tE_r(tE);
		  parlx[obsidx].set_piENE(piEN,piEE);
		  parlx[obsidx].fit_reinit();
		}
	    }


	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      fs[obsidx] = Event->fs[obsidx];
	    }

	  Event->deterror=0;

	  //Calculate the lightcurve
	  for(idx=0;idx<Event->nepochs;idx++)
	    {
	      obsidx = Event->obsidx[idx];
	      
	      tt = (Event->epoch[idx] - t0) / tE;
	      uu = u0;

	      if(pllx)
		{
		  tt += parlx[obsidx].tshift(idx-idxshift[obsidx]);
		  uu += parlx[obsidx].ushift(idx-idxshift[obsidx]);
		}

	      u = qAdd(tt,uu);
	      double amp2 = (u*u+2)/(u*sqrt(u*u+4)); 
	      if(u/rs>10&&u>5) amp = (u*u+2)/(u*sqrt(u*u+4)); 
	      else amp = wittFSMagnification(u,rs);
	      double amp3 = leeFSMagnification(u,rs);
	      //cerr << u << " " << amp << endl;
 
	      //store the results - remember F0=1 by definition
	      dF[pshift+idx] += dmult*fs[obsidx]*amp;

	      lc << param << " " << delta << " " << Event->epoch[idx] << " " << amp << " " << amp2 << " " << amp3 << " " << dmult*fs[obsidx]*amp << " " << dmult*fs[obsidx]*amp2 << " " << dmult*fs[obsidx]*amp << " " << tt << endl;
  
	    } //end for idx
	  lc << endl;
	} //end for delta

      //reset parallax parameters back to nominal values if necessary
      if(pllx && (param==3 || param==4 || param == 1))
	{
	  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
	    {
	      parlx[obsidx].set_tE_r(Event->tE_r);
	      parlx[obsidx].set_piENE(Event->piEN,Event->piEE);
	      parlx[obsidx].fit_reinit();
	    }
	}

    } //end for parameter

  //For the linear flux parameters we can work analytically
  //Do this numerically for the non-linear parameters
  int F0idx, fsidx;

  //Now calculate dF/dF0 and dF/dfs
  for(idx=0;idx<Event->nepochs;idx++)
    {
      obsidx = Event->obsidx[idx];
      F0idx = (Nparams+obsidx*2)*Event->nepochs+idx;
      fsidx = F0idx + Event->nepochs;
      dF[F0idx]=1 + fs[obsidx]*(Event->Atrue[idx]-1);
      dF[fsidx]=(Event->Atrue[idx]-1); //We have defined F0 as 1
    }     

  char fmfname[1000];

  if(Paramfile->choosefield<0)
    {
      sprintf(fmfname, "%s%s_%d_%d.det.fm", Paramfile->outputdir,
	      Paramfile->run_name, Event->instance, Event->id);
    }
  else
    {
      sprintf(fmfname, "%s%s_%d_%d_%d.det.fm", Paramfile->outputdir,
	      Paramfile->run_name, Event->instance, Paramfile->choosefield, 
	      Event->id);
    }


 

  
  fmfile_ptr = fopen(fmfname,"w");

  



  //Now calculate the bij matrix
  gsl_set_error_handler_off();
  gsl_matrix* bij = gsl_matrix_calloc(Ntotparams,Ntotparams); //inverse of cov
  gsl_matrix* cij = gsl_matrix_alloc(Ntotparams,Ntotparams); //cov
  gsl_permutation* permutation = gsl_permutation_alloc(Ntotparams);
  int signum;

  for(int i=0;i<Ntotparams;i++)
    {
      int ishift = i*Event->nepochs;
      for(int j=0;j<Ntotparams;j++)
	{
	  int jshift = j*Event->nepochs;

	  double sum=0;

	  //if(j<=i)
	  //  {
	  
	      for(idx=0;idx<Event->nepochs;idx++)
		{
		  sum += dF[ishift+idx] * dF[jshift+idx] / sqr(Event->Aerr[idx]);
		}
	      
	      gsl_matrix_set(bij,i,j,sum);
	      //  }
	      //else gsl_matrix_set(bij,i,j,gsl_matrix_get(bij,j,i));
	  
	}
    }

  fprintf(fmfile_ptr,"bij:\n");
  for(int i=0;i<Ntotparams;i++)
    {
      for(int j=0;j<Ntotparams;j++)
	{
	  fprintf(fmfile_ptr,"%6g\t",gsl_matrix_get(bij,i,j));
	}
      fprintf(fmfile_ptr,"\n");
    }
  fprintf(fmfile_ptr,"\n");

  //compute the LU decomposition
  gsl_linalg_LU_decomp(bij,permutation,&signum);

  fprintf(fmfile_ptr,"LU:\n");
  for(int i=0;i<Ntotparams;i++)
    {
      for(int j=0;j<Ntotparams;j++)
	{
	  fprintf(fmfile_ptr,"%6g\t",gsl_matrix_get(bij,i,j));
	}
      fprintf(fmfile_ptr,"\n");
    }
  fprintf(fmfile_ptr,"\n");

  //compute the inverse
  gsl_linalg_LU_invert(bij,permutation,cij);

  //print out the matrix
  fprintf(fmfile_ptr,"cij:\n");

   for(int i=0;i<Ntotparams;i++)
    {
      for(int j=0;j<Ntotparams;j++)
	{
	  fprintf(fmfile_ptr,"%6g\t",gsl_matrix_get(cij,i,j));
	}
      fprintf(fmfile_ptr,"\n");
    }
   fprintf(fmfile_ptr,"\n");


   fprintf(fmfile_ptr,"errors:\n");
   

   if(pllx)
     {
       fprintf(fmfile_ptr,"t0\tlogtE\tu0\tlogrs\tpiEN\tpiEE\t{F0\tlogfs}\n");
       fprintf(fmfile_ptr,"%6g\t%6g\t%6g\t%6g\t%6g\t%6g\t", Event->t0, 
	       logtE, Event->u0, log10(Event->rs), piEN,piEE);
     }
   else
     {
       fprintf(fmfile_ptr,"t0\tlogtE\tu0\tlogrs\tF0\tlogfs\n");
       fprintf(fmfile_ptr,"%6g\t%6g\t%6g\t%6g\t\t", Event->t0, 
	       logtE, Event->u0, log10(Event->rs));
     }

   for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
     {
       fprintf(fmfile_ptr,"%6g\t%6g\t",1.0,fs[obsidx]);
     }
   fprintf(fmfile_ptr,"\n");

   for(int i=0;i<Ntotparams;i++)
     {
       if(i==0 || i==2 || (pllx && (i==3 || i==4)) || (i-Nparams)%2==0)
	 fprintf(fmfile_ptr,"%6g\t",sqrt(gsl_matrix_get(cij,i,i)));
       else
	 fprintf(fmfile_ptr,"%6g\t",log(10)*sqrt(gsl_matrix_get(cij,i,i)));
     }
   fprintf(fmfile_ptr,"\n");

   fclose(fmfile_ptr);


  gsl_matrix_free(bij);
  gsl_matrix_free(cij);
  gsl_permutation_free(permutation);
  
  delete[] step;
  delete[] dF;
  delete[] fs;

}

