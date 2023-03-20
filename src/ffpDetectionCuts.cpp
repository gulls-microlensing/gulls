#include "lightcurveFitter.h"

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  double chi2point, chi2model;
  int nconsec=0;
  double lastt=1e30;

  int N3sig=0;

  int obsidx;


  Event->detected=0;

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  //Determine if an event is detected

  for(int obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      double t0_diff =1E50;
      double t0_diff_min = 1E50;
      int t0_diff_min_idx=0;
      Event->currentgroup=obsgroup;
      for(int grpidx=0;grpidx<int(Event->obsgroups[obsgroup].size());grpidx++)
	{
	  obsidx=Event->obsgroups[obsgroup][grpidx];
	  //  for(int obsidx=0; obsidx<Paramfile->numobservatories;obsidx++)
	  //    {
      
	  for(int idx=Event->nepochsvec[obsidx];idx<Event->nepochsvec[obsidx+1];idx++)
	    {
	      //Are the points detected at > 3 sigma? (must be consecutive here)
	      chi2point = sqr((Event->Aobs[idx]-1.0)/Event->Aerr[idx]);
	      chi2model = sqr((Event->Aobs[idx]-Event->Atrue[idx])/Event->Aerr[idx]);
	      
	      t0_diff = (Event->epoch[idx]-Event->t0);

	      if(abs(t0_diff)<abs(t0_diff_min) && grpidx==0)
		{
		  t0_diff_min=t0_diff;
		  t0_diff_min_idx = idx;
		} 
	       
	      if(Event->nosat[idx])
		{
		  //store the total chi^2 in the right place
		  Event->PSPL[obsgroup].chisq += chi2point-chi2model;
	      
		  if(chi2point>9)
		    {
		      if(Event->epoch[idx]-lastt>0 && Event->epoch[idx]-lastt<1)
			{
			  nconsec++;
			  if(nconsec>N3sig) N3sig=nconsec;
			}
		      else
			{
			  nconsec=1;
			}
		      //store chi_{3+}^2 in the observatory chi^2 vector
		      Event->PSPL[obsgroup].chisqvec[obsidx] += chi2point-chi2model;
		    }
		  else
		    {
		      nconsec=0;
		    }
		}

	      lastt=Event->epoch[idx];
	    }
	}
      

      //Doing scaling calculations if in the paramfile
      if(Paramfile->error_scaling)
	{
	  //Find uniform scaling factor that would have \delta\chi^2 > 300
	  Event->scale_factor_300 = sqrt(Event->PSPL[obsgroup].chisq/300);
	  //Find the nth point away from t0 such that we can find the scaling required to reduce n3sigma below threshold
	  int t0_idx;
	  //First we will do n3sig>=3
	  /*In this case
	    t0
	    |          |  |       |
	    t1         t2         t3
	    -->tmin-t0<0
	    -->thrid furthest point is -1 away from tmin_idx
	  */
	  if (t0_diff_min<0)
	    {
	      t0_idx = t0_diff_min_idx -1;
	    }
	  /*In this case
	    t0     
	    |       |  |          |
	    t1         t2         t3
	    -->tmin-t0>0
	    -->thrid furthest point is +1 away from tmin_idx
	  */
	  else
	    {
	      t0_idx = t0_diff_min_idx +1;
	    }
	  //Then for a SNR of 3
	  Event->scale_factor_n3sig3 = (Event->Atrue[t0_idx]-1)/(3*Event->Atrueerr[t0_idx]);
	  
	  //Next  we will do n3sig>=6
	  /*In this case
	    t0
	    |          |          |  |       |         |          |
	    t1         t2         t3         t4        t5         t6
	    -->tmin-t0<0
	    -->sixth furthest point is +3 away from tmin_idx
	  */
	  if (t0_diff_min<0)
	    {
	      t0_idx = t0_diff_min_idx +3;
	    }
	  /*In this case
	    t0
	    |          |          |       |  |         |          |
	    t1         t2         t3         t4        t5         t6
	    -->tmin-t0>0
	    -->thrid furthest point is +1 away from tmin_idx
	  */
	  else
	    {
	      t0_idx = to_diff_min_idx -3;
	    }
	  //Then for a SNR of 3
	  Event->scale_factor_n3sig6 = (Event->Atrue[t0_idx]-1)/(3*Event->Atrueerr[t0_idx]);
	}
      
      //Store the number of 3 sigma detections in the flatlc flag
      Event->flatlc[obsgroup] = N3sig;

      //Fit a PS lightcurve
      Event->deterror=0;

      if(Event->PSPL[obsgroup].chisq>80	 && Event->flatlc[obsgroup] >= 3)
	{
	  Event->detected=1;
	}

      Event->flatchi2[obsgroup] = Event->flatlc[obsgroup];
      if(Event->flatlc[obsgroup]>0) Event->flatlc[obsgroup]=0;
      else Event->flatlc[obsgroup]=1;
    }

}
