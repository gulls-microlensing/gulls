#include "lightcurveFitter.h"
#include "detectionCuts.h"
#include "fisher.h"
#include "fisherinversion.h"
#include "split.h"
#include <gsl/gsl_matrix.h>

#include<iostream>

#define DEBUGVAR 1

using namespace std;

void detectionCuts(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses)
{
  
  stringstream ss;
  string fmfname, header, values;

  int obsgroup, grpidx, obsidx;

  double pllx=Paramfile->pllxMultiplyer;

  string params_nopllx=string("t0 log_tE u0 alpha log_s log_q log_rho");
  string params_pllx=string(" piEN piEE");
  vector<string> params_obs;
  split(string("Fbase fs"),params_obs);

  gsl_matrix* cij; //covariance matrix

  int nfixpar = 7 + (pllx?2:0);


  vector<string> parstrings;
  if(pllx) 
	split(params_nopllx+params_pllx,parstrings);
  else
	split(params_nopllx,parstrings);

  //dont bother if there was an error generating the lightcurve
  if(Event->lcerror) return;

  //Determine if an event is detected 

  //Fit a PS lightcurve
  Event->deterror=0;
  Event->detected=0;

  for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
    {
      Event->currentgroup=obsgroup;
      if(DEBUGVAR) cout << "Entering lightcurveFitter" << endl;
      lightcurveFitter(Paramfile, World, Event);
      if(DEBUGVAR) cout << "Leaving lightcurveFitter" << endl;

      //reset parallax parameters back to nominal values if necessary
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
												 Event->piE, Event->thE);
		  Event->pllx[obsidx].compute_tushifts();
		}

	  int nobs = int(Event->obsgroups[obsgroup].size());	  
	  int nparams = nfixpar + 2*nobs;

	  //Generate the header information
	  for(int i=0;i<nparams;i++)
		{
		  ss.str(string(""));
		  ss << "ObsGroup_" << obsgroup << "_";
		  if(i<nfixpar)
			{
			  ss << parstrings[i] << "_err";
			}
		  else
			{
			  //for(grpidx=0;grpidx<nobs;grpidx++)
			  grpidx = int((i-nfixpar)/2);
			  obsidx=Event->obsgroups[obsgroup][grpidx];
			  
			  if((i-nfixpar)%2==0)
				{
				  ss << "Fbase" << obsidx << "_err";
				}
			  else
				{
				  ss << "fs" << obsidx << "_err";
					}
			}
		  Event->dataheader.push_back(ss.str());
		}
	  

      /* If finite source star effects important, fit FS lightcurve */
      if(abs(Event->umin) < 10*Event->rs && Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared)
		{
		  Event->flag_needFS[obsgroup]=1;
		  if(DEBUGVAR) cout << "Entering lightcurveFitter_FS" << endl;
		  lightcurveFitter_FS(Paramfile, World, Event); 
		  if(DEBUGVAR) cout << "Leaving lightcurveFitter_FS" << endl;
		  
		  //reset parallax parameters back to nominal values if necessary
		  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
			{
			  Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
													 Event->piE, Event->thE);
			  Event->pllx[obsidx].compute_tushifts();
			}

		}

      if(Event->deterror) 
		{
		  fprintf(stderr,"\nDiscarding event %d (Failed fit FS, errval=%d) \n", 
				  Event->id, Event->deterror);
		}
      else
		{
		  //did we detect it?
		  if((!Event->flag_needFS[obsgroup] && Event->PSPL[obsgroup].chisq>Paramfile->minChiSquared)
			 || (Event->flag_needFS[obsgroup] && Event->FSPL[obsgroup].chisq>Paramfile->minChiSquared))
			Event->detected=1;
		}
    }

  if(!Event->detected)
	{

	  for(obsgroup=0; obsgroup<int(Event->obsgroups.size()); obsgroup++)
		{
		  int nobs = int(Event->obsgroups[obsgroup].size());
		  int nparams = nfixpar + 2*nobs;
		  //Put dummy data in output
		  for(int i=0;i<nparams;i++)
			{
			  Event->data.push_back(9999);
			}
		}
	}
  else
    {
	  
      if(DEBUGVAR) cout << "Entering fisherMatrix" << endl;
      fisherMatrix(Paramfile, Event, World, Sources, Lenses);
      if(DEBUGVAR) cout << "Leaving fisherMatrix" << endl;

      //reset parallax parameters back to nominal values if necessary
      for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
		{
		  Event->pllx[obsidx].provide_murel_h_lb(Event->murel_l,Event->murel_b,
												 Event->piE, Event->thE);
		  Event->pllx[obsidx].compute_tushifts();
		}

      

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
		  ss << "t0\tlog(tE)\tu0\talpha\tlog(s)\tlog(q)\tlog(rho)\t";
		  if(pllx) ss << "piEN\tpiEE\t";
		  ss << "{F0\tfs}";
		  header = ss.str();
		  
		  ss.str(string(""));
		  ss.precision(12);
		  ss << scientific;
		  ss << Event->t0 << "\t" << log10(Event->tE_r) << "\t" 
			 << Event->u0 << "\t" << Event->alpha << "\t" 
			 << log10(Event->params[SS]) << "\t" 
			 << log10(Event->params[QQ]) << "\t" << log10(Event->rs) 
			 << "\t";
		  if(pllx) ss << Event->piEN << "\t" << Event->piEE << "\t";
		  for(grpidx=0;grpidx<nobs;grpidx++)
			{
			  obsidx = Event->obsgroups[obsgroup][grpidx];
			  ss << 1.0 << "\t" << Event->fs[obsidx] << "\t";
			} //for grpidx
		  values = ss.str();
	  
		  //Compile the differential vector
		  vector<double> dF;
		  vector<double> err;
	  
		  int nepochs=0;
		  for(grpidx=0;grpidx<nobs;grpidx++)
			{
			  obsidx=Event->obsgroups[obsgroup][grpidx];
			  int startidx = Event->nepochsvec[obsidx];
			  int endidx = Event->nepochsvec[obsidx+1];
			  nepochs += (endidx - startidx);
			} //for grpidx
	  
		  if(Paramfile->verbosity>1) cout << "\nFisher matrix, group " << obsgroup << endl;
		  
		  //int nobs = int(Event->obsgroups[obsgroup].size());
		  int nparams = nfixpar + 2*nobs;
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
				  //if(param==0) nepochs += endidx - startidx;
				  if(Paramfile->verbosity>2) cout << "nparams, nepochs, param, nobs, ofallparams, allshift, startidx, endidx, dFsize, Aerrsize = " << nparams << " " << nepochs << " " << param << " " << nobs << " " << ofallparams << " " << allshift << " " << startidx << " " << endidx << " " << Event->dF.size() << " " << Event->Aerr.size() << endl;
		  
				  for (int idx=startidx;idx<endidx;idx++)
					{
					  dF.push_back(Event->dF[idx + allshift]);
					  if(param==0) err.push_back(Event->Aerr[idx + allshift]);
					}
				  if(Paramfile->verbosity>2) cerr << "df+Aerr copied for parameter " << param << endl;
				} //for grpidx
			} //for param
	  
		  if(Paramfile->verbosity>1)
			{
			  cout << "Computing fisher inversion for group " << obsgroup << " with " << nparams << " params and " << nepochs << " datapoints" << endl << endl;
			}
	      

		  if(Paramfile->verbosity>2)
			{
			  cerr << "dF.size " << dF.size() << endl;
			  cerr << "err.size " << err.size() << endl;
			  cerr << "nparams " << nparams << endl;
			  cerr << "nepochs " << nepochs << endl;
			}

		  cij = gsl_matrix_alloc(nparams,nparams);
		  
		  fisherInversion(dF, err, nparams, nepochs, fmfname, header, values,
						  &cij);


		  //Put the data in output
		  for(int i=0;i<nparams;i++)
			{
			  Event->data.push_back(sqrt(gsl_matrix_get(cij,i,i)));
			}
		  
		  gsl_matrix_free(cij);
			   
		}
    }
}
