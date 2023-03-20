#include "lightcurveFitter.h"
#include "integerPowers.h"
#include "ephem.h"

#define DEBUGVAR 0

double my_f (const gsl_vector *v, void *params);
double fit_const(struct event *Event);

int lightcurveFitter(struct filekeywords* Paramfile, struct event *Event, int enablePllx)
{
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
  
  /* Starting point */

  Event->flag_needFS[obsgroup]=0; //Might not be true. Tested in main()
  Event->PSPL[obsgroup].rS = 0.0; //Just to set this value to zero
  Event->PSPL[obsgroup].ld1 = 0.0; //Just to set this value to zero
  Event->PSPL[obsgroup].chisq = 1e300;
  Event->PSPL[obsgroup].pllx.clear();
  for(int i=0;i<MAX_NUM_OBSERVATORIES;i++)
    {
      Event->PSPL[obsgroup].pllx.push_back(Event->pllx[i]);
      Event->PSPL[obsgroup].chisqvec[i]=0;
    }


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

  if(DEBUGVAR) printf("lightcurveFitter: nepochs = %d\n",Event->nepochs);

  if(Event->allsat)
    {
      Event->PSPL[obsgroup].chisq = 0;
      Event->flatchi2[obsgroup] = 0;
      delete[] fsblChi2obs;
      return GSL_SUCCESS;
    }

  Event->flatchi2[obsgroup] = fit_const(Event) - fsblChi2;

  if(Event->flatchi2[obsgroup]<60)
    {
      Event->flatlc[obsgroup]=1;
      Event->PSPL[obsgroup].chisq=0;
      delete[] fsblChi2obs;
      return GSL_SUCCESS;
    }
  else Event->flatlc[obsgroup]=0;

  if(Paramfile->pllxMultiplyer!=0&&enablePllx)
    {
      minex_func.n = 5;
      x = gsl_vector_alloc (5);
      ss = gsl_vector_alloc (5);

      gsl_vector_set (x, 3, Event->piEN);
      gsl_vector_set (x, 4, Event->piEE);

      gsl_vector_set (ss,3, 0.01);
      gsl_vector_set (ss,4, 0.01);
    }
  else
    {
      minex_func.n = 3;
      x = gsl_vector_alloc (3);
      ss = gsl_vector_alloc (3);
    }

  gsl_vector_set (x, 0, Event->u0);
  gsl_vector_set (x, 1, Event->t0);
  gsl_vector_set (x, 2, Event->tE_r);

  /* Set initial step sizes to 1 */
      
  gsl_vector_set (ss,0, 0.01);
  gsl_vector_set (ss,1, 10);
  gsl_vector_set (ss,2, 1);

    
  /* Initialize method and iterate */
  minex_func.f = &my_f;
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
	  if(Paramfile->pllxMultiplyer!=0&&enablePllx)
	    {
	      printf ("%5d %10.3e %10.3ef %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", 
		      (int)iter,
		      gsl_vector_get (s->x, 0), 
		      gsl_vector_get (s->x, 1), 
		      gsl_vector_get (s->x, 2), 
		      gsl_vector_get (s->x, 3), 
		      gsl_vector_get (s->x, 4), 
		      s->fval, size);
	    }
	  else
	    {
	      printf ("%5d %10.3e %10.3ef %10.3e f() = %7.3f size = %.3f\n",
		      (int)iter,
		      gsl_vector_get (s->x, 0), 
		      gsl_vector_get (s->x, 1), 
		      gsl_vector_get (s->x, 2), 
		      s->fval, size);
	    }
	}
    }
  while (status == GSL_CONTINUE && iter < 50000);
       
  gsl_vector_free(x);
  gsl_vector_free(ss);
  Event->PSPL[obsgroup].umin =      gsl_vector_get (s->x, 0);
  Event->PSPL[obsgroup].t0 =        gsl_vector_get (s->x, 1);
  Event->PSPL[obsgroup].tE =        gsl_vector_get (s->x, 2);
  if(Paramfile->pllxMultiplyer!=0&&enablePllx)
    {
      Event->PSPL[obsgroup].piEN =        gsl_vector_get (s->x, 3);
      Event->PSPL[obsgroup].piEE =        gsl_vector_get (s->x, 4);
    }
  Event->PSPL[obsgroup].chisq =     s->fval - fsblChi2;
  for(int i=0;i<MAX_NUM_OBSERVATORIES;i++) 
    Event->PSPL[obsgroup].chisqvec[i]-=fsblChi2obs[i];

  delete[] fsblChi2obs;
  gsl_multimin_fminimizer_free (s);

  if(DEBUGVAR) printf("lightcurve fitter: status = %d\n",status);

  return status;

}

double my_f(const gsl_vector *v, void *params)
{
  double sum(double x[],int m,int n);
  struct event *EventL = (struct event *)params;
  double umin, t0, tE, piEN, piEE, piE;
  double* Fu = new double[MAX_NUM_OBSERVATORIES];
  double* Fl = new double[MAX_NUM_OBSERVATORIES];
  double* chisqvec = new double[MAX_NUM_OBSERVATORIES];

  int obsgroup=EventL->currentgroup;
  int grpidx;

  double u2,det;       
  int idx,startidx,endidx,obsidx;
  double chisq = 0.0;
  double Chisq=0.0;

  double* A = new double[EventL->nepochs];
  double R1;
  double R2;

  double inverr_sqsum, R1sum, R1_sqsum, R2sum, R3sum, R4sum;

  double uu, tt;

  bool pllx=(v->size==3?false:true);

  for(obsidx=0;obsidx<MAX_NUM_OBSERVATORIES;obsidx++) chisqvec[obsidx]=0;
 
  umin = gsl_vector_get(v, 0);
  t0 = gsl_vector_get(v, 1);
  tE = gsl_vector_get(v, 2);
  piE=0; piEN=0; piEE=0;  //take care of the non-parallax case

  if(pllx)
    {
      piEN = gsl_vector_get(v, 3);
      piEE = gsl_vector_get(v, 4);
      piE = qAdd(piEE,piEN);
      
      for(grpidx=0;grpidx<int(EventL->obsgroups[obsgroup].size());grpidx++)
	{
	  obsidx=EventL->obsgroups[obsgroup][grpidx];
	  EventL->PSPL[obsgroup].pllx[obsidx].set_piENE(piEN,piEE);
	  EventL->PSPL[obsgroup].pllx[obsidx].set_tE_r(tE);
	  EventL->PSPL[obsgroup].pllx[obsidx].fit_reinit();
	}
    }

  double Aerr;
  double AA;
  double Aobs;
  double inverr2;

  //for(obsidx=0;obsidx<EventL->numobservatories;obsidx++)
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
		      //tt += EventL->PSPL[obsgroup].pllx[obsidx].tshift(idx-startidx);
		      //uu += EventL->PSPL[obsgroup].pllx[obsidx].ushift(idx-startidx);
		      tt += EventL->PSPL[obsgroup].pllx[obsidx].tshift(EventL->jdepoch[idx]);
		      uu += EventL->PSPL[obsgroup].pllx[obsidx].ushift(EventL->jdepoch[idx]);
		    }

		  u2 = sqr(tt) + sqr(uu);

		  A[idx] = AA = (u2 + 2 )/sqrt(u2*(u2+4.0));
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
	  for (idx=startidx;idx<endidx;idx++)
	    {
	      EventL->Afit[idx] = oFl*A[idx] + oFu;
	      if(EventL->nosat[idx])
		{
		  chisq += sqr( (EventL->Aobs[idx] - EventL->Afit[idx])
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

  if(Chisq < EventL->PSPL[obsgroup].chisq)
    {
      EventL->PSPL[obsgroup].chisq = Chisq;
      for(obsidx=0;obsidx<EventL->numobservatories;obsidx++)
	{
	  EventL->PSPL[obsgroup].chisqvec[obsidx] = chisqvec[obsidx];
	  EventL->PSPL[obsgroup].Fu[obsidx] = Fu[obsidx];
	  EventL->PSPL[obsgroup].Fl[obsidx] = Fl[obsidx];
	}
    }

  if(DEBUGVAR) printf("Chisq = %f\n",Chisq);

  //free up memory
  delete[] A;
  delete[] Fu;
  delete[] Fl;
  delete[] chisqvec;
 
  return Chisq;
 
}
     
double sum(double x[],int m,int n)
{
  double S= 0.0;
  for (int ind = m;ind<n;ind++) S+= x[ind];
  return(S);
}

double fit_const(struct event *Event)
{
  //fit a constant basline to the data to avoid fitting when events are missed
  //best fit is the weighted average

  double Asum;
  double Esum;
  double C;
  double invs2;
  double chi2=0;

  int idx,obsidx;
  
  int obsgroup = Event->currentgroup;
  int grpidx;

  for(grpidx=0;grpidx<int(Event->obsgroups[obsgroup].size());grpidx++)
    {
      obsidx=Event->obsgroups[obsgroup][grpidx];

      int startidx = Event->nepochsvec[obsidx];
      int endidx = Event->nepochsvec[obsidx+1];

      Asum=0;
      Esum=0;
	  
      //find the average
      if(startidx!=endidx && !Event->allsatobs[obsidx])
	{
	  for (idx=startidx;idx<endidx;idx++)
	    {
	      invs2 = 1.0/sqr(Event->Aerr[idx]);
	      Asum += Event->Aobs[idx]*invs2;
	      Esum += invs2;
	    }

	  C = Asum/Esum;

	  //calculate chi^2
	  for (idx=startidx;idx<endidx;idx++)
	    {
	      if(Event->nosat[idx]) 
		chi2 += sqr((Event->Aobs[idx]-C)/Event->Aerr[idx]);
	    }
	  
	}
    }

  return chi2;
}
