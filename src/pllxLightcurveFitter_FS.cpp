#include "lightcurveFitter.h"
#include "integerPowers.h"
#include "parallax.h"
#include "wittFSPL.h"

#include<iostream>

#define DEBUGVAR 0


int lightcurveFitter_FS(struct filekeywords* Paramfile, struct event *Event)
{
  double my_f_FS (const gsl_vector *v, void *params);

  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  
  size_t iter = 0;
  int status;
  double size;

  int obsgroup = Event->currentgroup;
  int grpidx;
  int obsidx;
  int idx;

  gsl_set_error_handler_off();

  //calculate the chi^2 of the true model
  double fsblChi2=0;
  double* fsblChi2obs = new double[MAX_NUM_OBSERVATORIES];
  for(int i=0;i<MAX_NUM_OBSERVATORIES;i++) fsblChi2obs[i]=0;

  for(grpidx=0;grpidx<int(Event->obsgroups[obsgroup].size());grpidx++)
    {
      obsidx=Event->obsgroups[obsgroup][grpidx];

      int startidx = Event->nepochsvec[obsidx];
      int endidx = Event->nepochsvec[obsidx+1];
	  
      if(startidx!=endidx && !Event->allsatobs[obsidx])
	{
	  for (idx=startidx;idx<endidx;idx++)
	    {
	      if(Event->nosat[idx])
		{
		  double chi2_ = sqr((Event->Aobs[idx]-Event->Atrue[idx])/Event->Aerr[idx]);
		  fsblChi2 += chi2_;
		  fsblChi2obs[Event->obsidx[idx]] += chi2_;
		}
	    }
	}
    }



  /* Starting point */

  Event->FSPL[obsgroup].chisq = 1e300;
  Event->FSPL[obsgroup].pllx.clear();
  for(int i=0;i<MAX_NUM_OBSERVATORIES;i++)
    {
      Event->FSPL[obsgroup].chisqvec[i]=0;
      Event->FSPL[obsgroup].pllx.push_back(Event->pllx[i]);
    }

  if(DEBUGVAR) printf("lightcurveFitter: nepochs = %d\n",Event->nepochs);
  

  if(Paramfile->pllxMultiplyer!=0)
    {
      minex_func.n = 6;
      x = gsl_vector_alloc (6);
      ss = gsl_vector_alloc (6);

      gsl_vector_set (x, 4, Event->piEN);     
      gsl_vector_set (x, 5, Event->piEE);     

      gsl_vector_set (ss,4, 0.01);
      gsl_vector_set (ss,5, 0.01);
    }
  else
    {
      minex_func.n = 4;
      x = gsl_vector_alloc (4);
      ss = gsl_vector_alloc (4);
    }
      
  gsl_vector_set (x, 0, Event->u0);
  gsl_vector_set (x, 1, Event->t0);
  gsl_vector_set (x, 2, Event->tE_r);
  gsl_vector_set (x, 3, Event->rs); 

  /* Set initial step sizes */
  gsl_vector_set (ss,0, 0.01);
  gsl_vector_set (ss,1, 10);
  gsl_vector_set (ss,2, 1);
  gsl_vector_set (ss,3, 1e-4);
   
  /* Initialize method and iterate */
  minex_func.f = &my_f_FS;
  minex_func.params = (void *)Event;
  
  s = gsl_multimin_fminimizer_alloc (T, minex_func.n);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status)
	{
	  printf("Problem with gsl_multimin_fminimizer_iterate. status = %d\n",status); 
	  break;
	}
     
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-3);
           
      if(DEBUGVAR)
	{
	  if(Paramfile->pllxMultiplyer!=0)
	    {
	      printf ("%5d %10.3e %10.3ef %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		      (int)iter,
		      gsl_vector_get (s->x, 0), 
		      gsl_vector_get (s->x, 1), 
		      gsl_vector_get (s->x, 2), 
		      gsl_vector_get (s->x, 3),
		      gsl_vector_get (s->x, 4),
		      gsl_vector_get (s->x, 5), 
		      s->fval, size);
	    }
	  else
	    {
	      printf ("%5d %10.3e %10.3ef %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		      (int)iter,
		      gsl_vector_get (s->x, 0), 
		      gsl_vector_get (s->x, 1), 
		      gsl_vector_get (s->x, 2), 
		      gsl_vector_get (s->x, 3),
		      s->fval, size);
	    }
	}
    }
  while (status == GSL_CONTINUE && iter < 50000);
       
  gsl_vector_free(x);
  gsl_vector_free(ss);
  Event->FSPL[obsgroup].umin = gsl_vector_get (s->x, 0);
  Event->FSPL[obsgroup].t0 =   gsl_vector_get (s->x, 1);
  Event->FSPL[obsgroup].tE =   gsl_vector_get (s->x, 2);
  Event->FSPL[obsgroup].rS =   abs(gsl_vector_get (s->x, 3)); //Allow it to move to -ve rs
  if(Paramfile->pllxMultiplyer!=0)
    {
      Event->FSPL[obsgroup].piEN = gsl_vector_get (s->x, 4);
      Event->FSPL[obsgroup].piEE = gsl_vector_get (s->x, 5);
    }
  Event->FSPL[obsgroup].chisq = s->fval - fsblChi2;
  for(int i=0;i<MAX_NUM_OBSERVATORIES;i++) 
    Event->FSPL[obsgroup].chisqvec[i]-=fsblChi2obs[i];

  gsl_multimin_fminimizer_free (s);
  Event->flag_needFS[obsgroup] = 1; 

  delete[] fsblChi2obs;

  char lcfname[1000];
  FILE * lcfile_ptr;
  double u0, tE, t0, rs;
  double u, mu, Amp;
  int jdx;

  if(DEBUGVAR)
    {
      /* Output the best lightcurve */

      sprintf(lcfname,"fslc%d.txt",Event->id);
      lcfile_ptr = fopen(lcfname,"w");
      if(lcfile_ptr == NULL)
	{
	  printf("Unable to open fslc file: %s",lcfname);
	  return status;
	}

      u0 = Event->FSPL[obsgroup].umin;
      t0 = Event->FSPL[obsgroup].t0;
      tE = Event->FSPL[obsgroup].tE;
      rs = Event->FSPL[obsgroup].rS;

      for(jdx=0;jdx<Event->nepochs;jdx++)
	{
	  //u = sqrt(u0*u0 + pow((Event->epoch[jdx] - t0)/tE , 2));
	  u = qAdd(u0, (Event->epoch[jdx] - t0)/tE);
	  muVisibility(&mu, rs, u, 0.0);

	  Amp = Event->FSPL[obsgroup].Fu[Event->obsidx[jdx]] 
	    + mu*Event->FSPL[obsgroup].Fl[Event->obsidx[jdx]];
	  fprintf(lcfile_ptr,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
		  Event->epoch[jdx],Amp,mu,u,u0,t0,tE,rs);
	}

      fclose(lcfile_ptr);
    }

  return status;

}

double my_f_FS (const gsl_vector *v, void *params)
{
  double sum(double x[],int m,int n);
  struct event *EventL = (struct event *)params;
  double umin, t0, tE;
  double mu,rs,ld1;
  double piEN, piEE, piE;
  double u,det;
  double uu, tt;
  int idx,obsidx,startidx,endidx;
  double chisq = 0.0;
  double Chisq=0.0;
  double* A = new double[EventL->nepochs];
  double R1, R2, inverr2;

  double* Fu = new double[MAX_NUM_OBSERVATORIES];
  double* Fl = new double[MAX_NUM_OBSERVATORIES];
  double* chisqvec = new double[MAX_NUM_OBSERVATORIES];

  for(obsidx=0;obsidx<MAX_NUM_OBSERVATORIES;obsidx++) chisqvec[obsidx]=0;

  int obsgroup=EventL->currentgroup;
  int grpidx;

  bool pllx=(v->size==4?false:true);  

  umin = gsl_vector_get(v, 0);
  t0 = gsl_vector_get(v, 1);
  tE = gsl_vector_get(v, 2);
  rs = abs(gsl_vector_get(v, 3));  //Allow it to move to -ve rs
  piE=0; piEE=0; piEN=0; //take care of the non-parallax case
  if(pllx)
    {
      piEN = gsl_vector_get(v, 4);
      piEE = gsl_vector_get(v, 5);
      piE = qAdd(piEE,piEN);

      for(grpidx=0;grpidx<int(EventL->obsgroups[obsgroup].size());grpidx++)
	{
	  obsidx=EventL->obsgroups[obsgroup][grpidx];
	  EventL->FSPL[obsgroup].pllx[obsidx].set_piENE(piEN,piEE);
	  EventL->FSPL[obsgroup].pllx[obsidx].set_tE_r(tE);
	  EventL->FSPL[obsgroup].pllx[obsidx].fit_reinit();
	}
    }

  ld1 = 0.0;  

  double R1sum, R1_sqsum, R2sum, R3sum, R4sum, inverr_sqsum; 
  double Aobs, Aerr, AA;

  for(grpidx=0;grpidx<int(EventL->obsgroups[obsgroup].size());grpidx++)
    {
      obsidx=EventL->obsgroups[obsgroup][grpidx];
      startidx = EventL->nepochsvec[obsidx];
      endidx = EventL->nepochsvec[obsidx+1];

      R1_sqsum=0;
      R2sum=0;
      R3sum=0;
      R4sum=0;
      inverr_sqsum=0;

      if(endidx>startidx && !EventL->allsatobs[obsidx])
	{
	  for(idx=startidx;idx<endidx;idx++)
	    {
	      if(EventL->nosat[idx])
		{
		  Aobs = EventL->Aobs[idx];
		  Aerr = EventL->Aerr[idx];

		  tt = (EventL->epoch[idx] - t0) / tE;
		  uu = umin;

		  if(pllx)
		    {
		      //tt += EventL->FSPL[obsgroup].pllx[obsidx].tshift(idx-startidx);
		      //uu += EventL->FSPL[obsgroup].pllx[obsidx].ushift(idx-startidx);
		      tt += EventL->FSPL[obsgroup].pllx[obsidx].tshift(EventL->jdepoch[idx]);
		      uu += EventL->FSPL[obsgroup].pllx[obsidx].ushift(EventL->jdepoch[idx]);
		    }

		  u = qAdd(uu,tt);

		  muVisibility(&mu,  rs, u, ld1);
		  A[idx] = AA = mu;

		  inverr_sqsum += (inverr2 = 1.0/sqr(Aerr));
		  R1sum += (R1 = AA/Aerr);
		  R1_sqsum += sqr(R1);
		  R2sum += (R2 = Aobs*inverr2);
		  R3sum += AA*inverr2;
		  R4sum += AA*R2;
		}
	    }

	  det = 1.0/(inverr_sqsum*R1_sqsum - R3sum*R3sum);
   
	  double oFu, oFl;
	  /* F_u */
	  Fu[obsidx] = oFu = (R1_sqsum*R2sum - R3sum*R4sum)*det; 
	  /* F_l */
	  Fl[obsidx] = oFl = (R4sum*inverr_sqsum - R3sum*R2sum)*det;  

	  chisq=0.0;   
	  
	  for(idx=startidx;idx<endidx;idx++)
	    {
	      EventL->Afit[idx] = oFl*A[idx] + oFu;
	      if(EventL->nosat[idx])
		{
		  //removed + oFu from inner bracket on line below
		  chisq += sqr((EventL->Aobs[idx] - EventL->Afit[idx])
			       / EventL->Aerr[idx]);
		} 
	    }
   
	  chisqvec[obsidx] = chisq;
	  Chisq += chisq;
	}
      else
	{
	  Fu[obsidx] = -9999;
	  Fl[obsidx] = -9999;
	  chisqvec[obsidx]=0;
	}
	  
    }

  if(Chisq < EventL->FSPL[obsgroup].chisq)
    {
      EventL->FSPL[obsgroup].chisq = Chisq;
      for(obsidx=0;obsidx<EventL->numobservatories;obsidx++)
	{
	  EventL->FSPL[obsgroup].chisqvec[obsidx] = chisqvec[obsidx];
	  EventL->FSPL[obsgroup].Fu[obsidx] = Fu[obsidx];
	  EventL->FSPL[obsgroup].Fl[obsidx] = Fl[obsidx];
	}
    }

  //free up memory
  delete[] A;
  delete[] Fu;
  delete[] Fl;
  delete[] chisqvec;
 
  return Chisq;
 
}
     

void muVisibility(double* mu, double rs, double z0, double ld1)
{
  ld1=0; /*Don't include limb darkening yet*/

  *mu = wittFSMagnification(z0,rs);
  return;
}
