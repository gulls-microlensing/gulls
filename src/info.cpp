#include "info.h"
#include "constdefs.h"

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

#include "split.h"

#define DEBUGVAR 0

void printEpochs(struct obsfilekeywords World[], int obsidx)
{
  FILE *outfile_ptr1;
  char filename[1000];
  
  sprintf(filename,"epochs_%s",World[obsidx].name);

  outfile_ptr1 = fopen(filename,"w");
  int ind;
  char str[1100];
  if (outfile_ptr1 == NULL)
    {
      sprintf(str,"Unable to open output file: %s",filename);
      fmtline(str,WIDTH,"(info)");
      exit(1);
    }
  else sprintf(str,"Output file: %s",filename);
  fmtline(str,WIDTH,"(info)");

  for(ind=0;ind<World[obsidx].nepochs;ind++)
    {
      fprintf(outfile_ptr1,"%f %d %d\n",
	      World[obsidx].epoch[ind],World[obsidx].field[ind],
	      World[obsidx].weather[ind]);
    }

  fclose(outfile_ptr1);

}

void showIsseenby(struct obsfilekeywords World[], struct event *Event, int nobs)
{

  int ind,obsidx;
  char str[100];
  for(obsidx = 0;obsidx<nobs;obsidx++)
    {

      for(ind = 0;ind<World[obsidx].nfields;ind++)
	{
	  sprintf(str,"Observatory: %d Field: %d Seen: %d",
		  obsidx,ind,Event->isseenby[obsidx][ind]);
	  fmtline(str,45,"(showIsseenby)");
	}
    }
}

void writeHeader(struct filekeywords* Paramfile, struct event *Event, struct slcat* Sources, struct slcat* Lenses, ofstream& ofile)
{

  char paramstr[100];
  //position and event data - +6+1 = 1
  ofile << "EventID" << " " << "SubRun" << " ";
  ofile << "Field" << " ";
  ofile << "galactic_l" << " " << "galactic_b" << " " << "ra_deg" << " " 
		<< "dec_deg" << " ";
  //ofile << "| ";

  //source data - +6+1 = 8
  ofile << "SourceID" << " ";
  for(int i=0;i<Sources->datakey.size();i++)
	{
	  ofile << "Source_" << Sources->datakey[i] << " ";
	}
  //for(int i=0;i<sOutputCols;i++)
  //  {
  //    ofile << Sources->data[Event->source][sOutputColumns[i]] << " ";
  //  }

  //lens data - +10+1 = 15
  ofile << "LensID" << " ";
  for(int i=0;i<Lenses->datakey.size();i++)
	{
	  ofile << "Lens_" << Lenses->datakey[i] << " ";
	}
  //for(int i=0;i<lOutputCols;i++)
  //  {
  //    ofile << Lenses->data[Event->lens][lOutputColumns[i]] << " ";
  //  }
  //ofile << "| ";

  
  //microlensing paramters - 11+1 = 26
  ofile << "u0lens1" << " " << "alpha" << " ";
  ofile << "t0lens1" << " ";
  ofile << "tref" << " ";
  ofile << "tcroin" << " ";
  ofile << "ucroin" << " " << "rcroin" << " ";
  ofile << "tE_ref" << " " << "tE_helio" << " " << "rE" << " " << "thetaE" << " " 
		<< "rho" << " " << "piE" << " " << "piEN" << " " << "piEE" << " "
		<< "murel_helio_alpha" << " " << "murel_helio_delta" << " " << "murel_helio_l" << " " << "murel_helio_b" << " " << "murel_helio_lambda" << " " << "murel_helio_beta" << " " << "murel_helio" << " "
		<< "murel_ref_alpha" << " " << "murel_ref_delta" << " " << "murel_ref_l" << " " << "murel_ref_b" << " " << "murel_ref_lambda" << " " << "murel_ref_beta" << " " << "murel_ref" << " "
		<< "vtilde_helio" << " " << "vtilde_ref" << " " << "v_ref" << " "
		<< "vtilde_helio_N" << " " << "vtilde_helio_E" << " " <<  "vtilde_ref_N" << " " << "vtilde_ref_E" << " " << "v_ref_N" << " " << "v_ref_E" << " "
		<< "piEll" << " " << "piErp" << " "
		<< "vt" << " " << "LDgamma" << " ";
  //ofile << " | "; 

  //planet data - 7+1 = 38
  for(int i=0;i<NPLANETINPUT+NPLANETDERIV;i++)
    {
	  sprintf(paramstr,"%d",i);
	  if(Event->paramsHeader[i].size()>0) strcpy(paramstr,Event->paramsHeader[i].c_str());
      ofile << "Planet_" << paramstr << " ";
    }
  //ofile << "| ";
  
  //weights - 3+1 = 46
  ofile << "u0max" << " " << "t0range" << " " << "weight_scale" << " " << "raw_weight" << " " << "weight" << " ";
  //ofile << " | ";

  //simulation results

  //magnitudes - 10+2 = 50
  for(int i=0;i<Paramfile->Nfilters;i++)
    ofile << "Source_" << Sources->magkey[i] << " ";
  //ofile << "| ";
  for(int i=0;i<Paramfile->Nfilters;i++) 
    ofile << "Lens_" << Lenses->magkey[i] << " ";
  //ofile << "| ";

  //blending - 4+1 = 62
  for(int i=0;i<Paramfile->numobservatories;i++)
    ofile << "Obs_" << i << "_fs" << " ";
  //ofile << "| ";

  //simulation details
  ofile << "NumObsGroups" << " " << "ErrorFlag" << " ";
  //ofile << " | ";


  //chi^2 and obsgroup results
  for(int obsgroup=0;obsgroup<int(Event->obsgroups.size());obsgroup++)
    {
      ofile << "ObsGroup_" << obsgroup << "_" << "flatlc" << " ";
	  ofile << "ObsGroup_" << obsgroup << "_" << "flatchi2" << " ";
	  ofile << "ObsGroup_" << obsgroup << "_" << "FiniteSourceflag" << " ";
	  ofile << "ObsGroup_" << obsgroup << "_" << "chi2" << " ";
	  ofile << Event->obsgroupoutputheader[obsgroup] << " ";
    }
  

  // = 78

  //Generic data - determined by the user
  if(Event->data.size()>0)
    {
	  if(Event->dataheader.size()==Event->data.size())
		{
		  //ofile << " || ";
		  for(int i=0;i<int(Event->dataheader.size());i++)
			{
			  ofile << Event->dataheader[i] << " ";
			}
		}
	  else
		{
		  for(int i=0;i<int(Event->data.size());i++)
			{
			  ofile << "EventData" << i << " ";
			}
		}
    }

  if(Paramfile->error_scaling)
    {
	  //ofile << "| ";
	  ofile << "scale_factor_300" << " " << "scale_factor_n3sig3" << " " << "scale_factor_n3sig6" << " ";
    }

  //whether a lightcurve output file was generated
  ofile << "LCOutput" << " ";

  ofile << endl;

}

void writeEventParams(struct filekeywords* Paramfile, struct event *Event, struct slcat* Sources, struct slcat* Lenses, ofstream& ofile)
{
  ofile.precision(12);
  ofile << scientific;
  //position and event data - +6+1 = 1
  ofile << Event->id << " " << Event->instance << " ";
  ofile << Event->field << " ";
  ofile << Event->l << " " << Event->b << " " << Event->ra*TO_DEG << " " 
	<< Event->dec*TO_DEG << " ";

  //source data - +6+1 = 8
  int sn = Event->source;
  ofile << sn << " ";
  for(int i=0;i<Sources->data[sn].size();i++)
	{
	  ofile << Sources->data[sn][i] << " ";
	}
  //for(int i=0;i<sOutputCols;i++)
  //  {
  //    ofile << Sources->data[Event->source][sOutputColumns[i]] << " ";
  //  }
  //ofile << "| ";

  //lens data - +10+1 = 15
  int ln = Event->lens;
  ofile << ln << " ";
  for(int i=0;i<Lenses->data[ln].size();i++)
	{
	  ofile << Lenses->data[ln][i] << " ";
	}
  //for(int i=0;i<lOutputCols;i++)
  //  {
  //    ofile << Lenses->data[Event->lens][lOutputColumns[i]] << " ";
  //  }
  //ofile << "| ";

  //microlensing paramters - 11+1 = 26
  ofile << Event->u0 << " " << Event->alpha << " ";
  ofile << Event->t0 << " ";
  ofile << Paramfile->tref << " ";
  ofile << Event->tcroin << " ";
  ofile << Event->ucroin << " " << Event->rcroin << " ";
  double murel_l = Lenses->data[ln][MUL]-Sources->data[sn][MUL];
  double murel_b = Lenses->data[ln][MUB]-Sources->data[sn][MUB];
  double murel_ref = Event->thE/Event->tE_r * DAYINYR;


  ofile << Event->tE_r << " " << Event->tE_h << " " << Event->rE << " " << Event->thE << " " 
		<< Event->rs << " " << Event->piE << " " << Event->piEN << " " << Event->piEE << " "
		<< Event->murel*Event->pllx[0].mua_h << " " << Event->murel*Event->pllx[0].mud_h << " " << Event->murel*Event->pllx[0].mul_h << " " << Event->murel*Event->pllx[0].mub_h << " " << Event->murel*Event->pllx[0].mulam_h << " " << Event->murel*Event->pllx[0].mubet_h << " " << Event->pllx[0].murel_h << " "
		<< murel_ref*Event->pllx[0].mua_r << " " <<  murel_ref*Event->pllx[0].mud_r << " " << murel_ref*Event->pllx[0].mul_r << " " << murel_ref*Event->pllx[0].mub_r << " " << murel_ref*Event->pllx[0].mulam_r << " " << murel_ref*Event->pllx[0].mubet_r << " " << murel_ref*Event->pllx[0].murel_r << " "
		<< Event->pllx[0].vtilde_h << " " << Event->pllx[0].vtilde_r << " " << Event->pllx[0].v_ref << " "
		<< Event->pllx[0].vtilde_N_h << " " << Event->pllx[0].vtilde_E_h << " " <<  Event->pllx[0].vtilde_N_r << " " << Event->pllx[0].vtilde_E_r << " " << Event->pllx[0].v_ref_N << " " << Event->pllx[0].v_ref_E << " "
		<< Event->pllx[0].piEll << " " << Event->pllx[0].piErp << " "
		<< Event->vt << " " << Event->gamma << " "; 

  
  //planet data - 7+1 = 38
  for(int i=0;i<NPLANETINPUT+NPLANETDERIV;i++)
    {
      ofile << Event->params[i] << " ";
    }
  //ofile << "| ";
  
  //weights - 3+1 = 46
  ofile << Event->u0max << " " << Event->t0range << " " << Event->weight_scale << " " << Event->raww << " " << Event->w << "  ";

  //simulation results

  //magnitudes - 10+2 = 50
  for(int i=0;i<Paramfile->Nfilters;i++)
    ofile << Sources->mags[Event->source][i] << " ";
  //ofile << "| ";
  for(int i=0;i<Paramfile->Nfilters;i++) 
    ofile << Lenses->mags[Event->lens][i] << " ";
  //ofile << "| ";

  //blending - 4+1 = 62
  for(int i=0;i<Paramfile->numobservatories;i++)
    ofile << Event->fs[i] << " ";
  //ofile << "| ";

  //simulation details
  ofile << int(Event->obsgroups.size()) << " " << Event->allsat + 2*(!Event->nepochs) << " ";


  //chi^2 results
  for(int obsgroup=0;obsgroup<int(Event->obsgroups.size());obsgroup++)
    {
      //what is being done here? 
      ofile << Event->flatlc[obsgroup] << " " << Event->flatchi2[obsgroup] << " " << Event->flag_needFS[obsgroup] << " " << 
	(!Event->flag_needFS[obsgroup]?Event->PSPL[obsgroup].chisq:Event->FSPL[obsgroup].chisq) << " ";
      ofile << Event->obsgroupoutput[obsgroup] << " ";
    }
  
  // = 78

  //Generic data - determined by the user
  if(Event->data.size()>0)
    {
      for(int i=0;i<int(Event->data.size());i++)
		{
		  ofile << Event->data[i] << " ";
		}
    }
  

  if(Paramfile->error_scaling)
    {
        ofile  << Event->scale_factor_300 << " " << Event->scale_factor_n3sig3 << " " << Event->scale_factor_n3sig6 << " ";
    }

  //whether a lightcurve output file was generated
  ofile << Event->outputthis*Event->detected << " ";

  ofile << endl;
  
  
}

void writeEventLC(struct event *Event , int idx, char eventprefix[])
{

  FILE *lc_ptr;
  char outfile[1000];
  char str[1100];
  
  int jdx;
  
  sprintf(outfile,"%s.%06d",eventprefix,idx);

  lc_ptr = fopen(outfile,"w");
 
  if (lc_ptr == NULL)
    {
      sprintf(str,"Unable to open output file: %s",outfile);
      fmtline(str,2*WIDTH,"FAILED (info:writeEventLC)");
      exit(1);
    }

 for(jdx=0;jdx<Event->nepochs;jdx++)
   {
     if(Event->nosat[jdx])
       {
	 fprintf(lc_ptr,"%f %f %f %d\n",
		 Event->epoch[jdx],Event->Aobs[jdx],Event->Aerr[jdx],
		 Event->obsidx[jdx]);
       }
   }
 
 fclose(lc_ptr);

}

void errorHandler(int errval)
{

  /*   Error classes */

  /*   1xxx : In main code, in functions and subfunctions outside event
       generation loop */
  /*   2xxx : In main code, in functions and subfunctions inside loop */
  /*   3xxx : In mag.f and subfunctions  */

  char str[1000];

  switch (errval)
    {

    case 3001:
      sprintf(str,"limb darkening profile not yet implemented (mag.f)");
      fmtline(str,WIDTH,"(errorHandler)"); 
      sprintf(str,"Exiting");
      fmtline(str,WIDTH,"(errorHandler)"); 
      exit(1);
      break;

    case 3111:
      sprintf(str,
	      "coordinate system incorrect, m1 z2+m2 z1 !=0 (findImages.f)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      sprintf(str,"Exiting");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3112:
      sprintf(str,"more than 2 images already found should not enter findCloseImages (findImages.f:findCloseImages)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      sprintf(str,"Exiting");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3113:
      sprintf(str,
	      " image number incorrect in newImage (findImages.f:newImage)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      sprintf(str,"Exiting");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3114:
      sprintf(str,"nimages incorrect (findImages.f:checkImageProperties)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 

      break;
      
    case 3115:
      sprintf(str, "total magnification < 1, impossible (findImages.f:checkImageProperties)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 

      break;

    case 3201:
      sprintf(str,
	      " linear interpolation failed (uniform.f:uniformInterpolate)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3202:
      sprintf(str," failed in uniformInterpolate (uniform.f:uniformFunc)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;


    case 3211:
      sprintf(str," too many points in i & j (track.f:connectEnds)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;
      
    case 3212:
      sprintf(str," too many points in adaptive (track.f:adaptiveFindcomplextracks)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;
      
    case 3213:
      sprintf(str," too many points in i & j (track.f:jumpConnects)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;
    
    case 3214:
      sprintf(str," complex root finding failed (track.f:complexSolve)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3215:
      sprintf(str," wrong image number (track.f:complexSolve)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 3216:
      sprintf(str," too many points+1 (track.f:connectSegments)");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    case 4001:
      sprintf(str," error from backup generator ");
      fmtline(str,2*WIDTH,"(errorHandler)"); 
      break;

    default:
      sprintf(str,"Unknown error code: %d (info.c:errorHandler)",errval);
      fmtline(str,WIDTH,"FATAL"); 
      exit(1);
    }

}


void progressbar(int idx, int niter, time_t st1)
{
  char bar[1000];
  sprintf(bar,"\r[");
  int nchars,ind;
  char tmp[1000];
  double frac;
  double etc;
  frac = (double)(idx+1)/(double)niter;
  int hrs,mins;
  float sec;

  etc=difftime(time(0),st1) * (1/frac -1);
  
  nchars = int(floor(frac*40.0));
  for(ind=0;ind<nchars;ind++)
    {
      strcat(bar,"=");
    }
  for(ind=nchars;ind<40;ind++)
    {
      strcat(bar," ");
    }
  
  hrs = int(floor(etc/3600));
  mins = int(floor((etc - hrs*3600)/60));
  sec = etc - mins*60 - hrs*3600;
  
  sprintf(tmp,"] (%d/%d) ETC:%3dh %2dm %2ds",idx+1,niter,hrs,mins,(int)sec);
  
  strcat(bar,tmp);
  printf("%s",bar);

}

void clock2str(time_t now, time_t start)
{
  int hrs,mins;
  double sec, etc;
  etc = difftime(now,start);

  hrs = int(floor(etc/3600));
  mins = int(floor((etc - hrs*3600)/60));
  sec = etc - mins*60 - hrs*3600;
  
  printf("Execution time:  %3dh %2dm %2ds\n",hrs,mins,(int)sec);
  
}
