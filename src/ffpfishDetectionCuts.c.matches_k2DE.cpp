#include "lightcurveFitter.h"
#include "fisher.h"
#include "fisherinversion.h"
#include "mderrors.h"

#include<iostream>
#include<sstream>

#define DEBUGVAR 0

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{

  int obsgroup, grpidx, obsidx;
  double pllx=Paramfile->pllxMultiplyer;

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  double chi2point, chi2model;

  stringstream ss;
  string fmfname;
  string header;
  string values;

  Event->deterror=0;
  Event->detected=0;
  
  //Determine if an event is detected

  for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      Event->currentgroup=obsgroup;
      Event->obsgroupoutput[obsgroup] = string("999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 999999 ");

      Event->flag_needFS[obsgroup]=0;
      Event->PSPL[obsgroup].chisq=0;

      int N3sig=0;
      int nconsec=0;
      double lastt=1e30;

      for(grpidx=0;grpidx<int(Event->obsgroups[obsgroup].size());grpidx++)
	{
	  obsidx=Event->obsgroups[obsgroup][grpidx];
	  
	  for(int idx=Event->nepochsvec[obsidx];idx<Event->nepochsvec[obsidx+1];idx++)
	    {
	      //Are the points detected at > 3 sigma? (must be consecutive here)
	      chi2point = sqr((Event->Aobs[idx]-1.0)/Event->Aerr[idx]);
	      chi2model = sqr((Event->Aobs[idx]-Event->Atrue[idx])/Event->Aerr[idx]);

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

      //Store the number of 3 sigma detections in the flatlc flag
      Event->flatlc[obsgroup] = N3sig;

      if(Event->PSPL[obsgroup].chisq>80 && Event->flatlc[obsgroup] >= 3)
	{
	  Event->detected=1;
	}

      Event->flatchi2[obsgroup] = Event->flatlc[obsgroup];
      if(Event->flatlc[obsgroup]>0) Event->flatlc[obsgroup]=0;
      else Event->flatlc[obsgroup]=1;

    }

  if(Event->detected)
    {
      if(DEBUGVAR) cout << "Entering fisherMatrix" << endl;
      fisherMatrix(Paramfile, Event, World, Sources, Lenses);
      if(DEBUGVAR) cout << "Leaving fisherMatrix" << endl;

      for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
	{
	  Event->currentgroup=obsgroup;
	  int nobs = int(Event->obsgroups[obsgroup].size());

	  //Generate the header information
	  ss.str(string("")); 
	  ss << Paramfile->outputdir << Paramfile->run_name << "_" 
	     << Event->instance << "_";
	  if(Paramfile->choosefield>=0) ss << Paramfile->choosefield << "_";
	  ss << Event->id << ".det.fm." << obsgroup;
	  fmfname = ss.str();

	  if(!Event->outputthis) fmfname=string("");

	  ss.str(string("")); 
	  ss << "t0\ttE\tu0\trs\t";
	  if(pllx) ss << "piEN\tpiEE\t";
	  ss << "{F0\tfs}";
	  header = ss.str();
	  
	  ss.str(string(""));
	  //ss.precision(16);
	  ss << Event->t0 << "\t" << Event->tE_r  << "\t" << Event->u0 << "\t" 
	     << Event->rs << "\t";
	  if(pllx) ss << Event->piEN << "\t" << Event->piEE << "\t";	  
	  for(grpidx=0;grpidx<nobs;grpidx++)
	    {
	      obsidx = Event->obsgroups[obsgroup][grpidx];
	      ss << 1.0 << "\t" << Event->fs[obsidx] << "\t";
	    }
	  values = ss.str();

	  //Compile the differential vector
	  vector<double> dF;
	  vector<double> err;

	  int nfixpar = 4 + (pllx?2:0);
	  int nparams = nfixpar + 2*nobs;

	  int nepochs=0;
	  for(grpidx=0;grpidx<nobs;grpidx++)
	    {
	      obsidx=Event->obsgroups[obsgroup][grpidx];
	      int startidx = Event->nepochsvec[obsidx];
	      int endidx = Event->nepochsvec[obsidx+1];
	      nepochs += (endidx - startidx);
	    }

	  if(Paramfile->verbosity>1) cout << "\nFisher matrix, group " << obsgroup << endl;

	  for(int param=0;param<nparams;param++)
	    {
	      int parobsidx = Event->obsgroups[obsgroup][(param-nfixpar)/2];
	      int ofallparams = (param<nfixpar ? 
				 param : 
				 nfixpar + 2*parobsidx + (param-nfixpar)%2);
	      int allshift = Event->nepochs*ofallparams;

	      for(grpidx=0;grpidx<nobs;grpidx++)
		{
		  obsidx=Event->obsgroups[obsgroup][grpidx];
		  
		  int startidx = Event->nepochsvec[obsidx];
		  int endidx = Event->nepochsvec[obsidx+1];

		  if(Paramfile->verbosity>2) cout << "nparams, nepochs, param, nobs, ofallparams, allshift, startidx, endidx, dFsize, Aerrsize = " << nparams << " " << nepochs << " " << param << " " << nobs << " " << ofallparams << " " << allshift << " " << startidx << " " << endidx << " " << Event->dF.size() << " " << Event->Aerr.size() << endl;

		  for (int idx=startidx;idx<endidx;idx++)
		    {
		      //if(idx+allshift>=Event->dF.size()) cout << "Here's the segfault" << endl;
		      dF.push_back(Event->dF[idx + allshift]);
		      if(param==0) err.push_back(Event->Aerr[idx]);
		    } //for idx
		  if(Paramfile->verbosity>2) cout << "dF+Aerr copied for parameter " << param << endl;

		} //for grpidx
	    } //for param

	  if(Paramfile->verbosity>1)
	    {
	      cout << "Computing fisher inversion for group " << obsgroup << " with " << nparams << " params and " << nepochs << " datapoints" << endl;
	    }

	  gsl_matrix* cov = gsl_matrix_alloc(1,1);
	  fisherInversion(dF, err, nparams, nepochs, fmfname, header, values, 
			  &cov);

	  if(Paramfile->verbosity>1)
	    {
	      cout << "Transforming parameters to mass/pirel"  << endl << endl;
	    }

	  vector<double> bestmdcov(7,999999);
	  vector<int> bestfilters(2,999999);
	  vector<double> piEdir(3,999999);

	  if(Event->obsgroups[obsgroup].size()>=2)
	    {
	      piEdir = piEdirerrors(cov, nparams, Event, Sources, Lenses);

	      for(int obs1=0;obs1<int(Event->obsgroups[obsgroup].size());obs1++)
		{
		  for(int obs2=0;obs2<int(Event->obsgroups[obsgroup].size());obs2++)
		    {
		      int filter1 = World[Event->obsgroups[obsgroup][obs1]].filter;
		      int filter2 = World[Event->obsgroups[obsgroup][obs2]].filter;
		      int worstalpha=0;
		      double k2alpha=10; //only use approved filters!
		      int goodfilter=0;
		      
		      //alphas from fitting ax+b to Dartmouth fehp00afep2 
		      //isochrone for logg<4.5. Assumes that redder filter
		      //sets the zero magnitude angular diameter
		      if(filter1==1 && filter2==4)
			{
			  //g-z
			  k2alpha = 0.217;
			  goodfilter=0;
			}
		      if(filter1==1 && filter2==3)
			{
			  //g-i
			  k2alpha = 0.294;
			  goodfilter=0;
			}
		      if(filter1==1 && filter2==2)
			{
			  //g-r
			  k2alpha = 0.511;
			  goodfilter=0;
			}
		      if(filter1==2 && filter2==3)
			{
			  //r-i
			  k2alpha = 0.962;
			  goodfilter=1;
			}
		      if(filter1==2 && filter2==4)
			{
			  //r-z
			  k2alpha = 0.523;
			  goodfilter=1;
			}
		      if(filter1==3 && filter2==4)
			{
			  //i-z
			  k2alpha = 1.375;
			  goodfilter=1;
			}
		      if(filter1==6 && filter2==8) 
			{
			  //V-I
			  k2alpha = 0.427; 
			  goodfilter=1;
			}
		      

		      //vector<double> mdcov = mderrors(fmfname, Event->obsgroups[obsgroup][obs1], Event->obsgroups[obsgroup][obs2], k2alpha, 1, Event, Sources, Lenses);
		      vector<double> mdcov;
		      if(goodfilter)
			{
			  mdcov = mderrors(cov, nparams, Event->obsgroups[obsgroup][obs1], Event->obsgroups[obsgroup][obs2], k2alpha, 1, Event, Sources, Lenses);

			  if(mdcov[0]<bestmdcov[0] && mdcov[0]>0 && mdcov[2]>0)
			    {
			      bestmdcov = mdcov;
			      bestfilters[0] = Event->obsgroups[obsgroup][obs1];
			      bestfilters[1] = Event->obsgroups[obsgroup][obs2];
			    } //if better filter pair

			} //if goodfilter
		    } //for obs2
		} //for obs1

	      //cout << "mderr group " << obsgroup << " " << bestmdcov[0] << " " << bestmdcov[1] << " " << bestmdcov[2] << " " << bestfilters[0] << " " << bestfilters[1] << endl;
	    } //if 2 obs

	  stringstream mdss;
	  mdss << bestmdcov[0] << " " << bestmdcov[1] << " " << bestmdcov[2]  << " " << bestfilters[0] << " " << bestfilters[1] << " " << bestmdcov[3] << " " << bestmdcov[4] << " " << bestmdcov[5] << " " << bestmdcov[6] << " " << gsl_matrix_get(cov,3,3)/Event->rs/ln10 << " " << gsl_matrix_get(cov,4,4)/abs(Event->piEN)/ln10 << " " << gsl_matrix_get(cov,5,5)/abs(Event->piEE)/ln10 << " ";
	  Event->obsgroupoutput[obsgroup] = mdss.str();

	  //cout << "before free" << endl;
	  gsl_matrix_free(cov);
	  //cout << "after free" << endl;

	} //for obsgroup

	  

    } //if detected
  
}
