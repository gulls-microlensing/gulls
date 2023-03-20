#include "outputLightcurve.h"
#include "lightcurveFitter.h"
#include "zodiacalLight.h"
#include "astroFns.h"

#define DEBUGVAR 1

void outputLightcurve(struct event *Event, struct filekeywords* Paramfile, struct slcat* Sources, struct slcat* Lenses)
{
  void muVisibility(double *mu, double rs, double z0, double ld1);
  char lcfname[1000];
  FILE *lcfile_ptr;
  int fileOpen=0;
  int i, obsidx;
  char data[1000];
  char tmp[30];
  char extension[4];

  double t;
  
  //double u0_ps,t0_ps,tE_ps,t2_ps,u_ps,psmag,psamp;
  //double u0_fs,t0_fs,tE_fs,t2_fs,u_fs,fsmag,fsamp;

  if(!Paramfile->outputLightcurve 
     || !(Paramfile->outputOnErr || Paramfile->outputOnDet 
	  || Paramfile->outputOnAll)) return;
  else
    {
      if((Event->lcerror+Event->deterror))
	{
	  if(!Paramfile->outputOnErr) return;
	}
      else if((Event->detected))
	{
	  if(!Paramfile->outputOnDet) return;
	}
      else if((!Event->detected))
	{
	  if(!Paramfile->outputOnAll) return;
	}
    }

  if(Event->detected) strcpy(extension,"det");
  else if(Event->lcerror||Event->deterror) strcpy(extension,"err");
  else strcpy(extension,"all");

  sprintf(lcfname, "%s%s_%d_%d.%s.lc", Paramfile->outputdir,
	  Paramfile->run_name, Event->instance, Event->id, extension);

  if(DEBUGVAR) printf("lcname: %s\n",lcfname);
  
  if(Event->nepochs>0)
    {
      lcfile_ptr = fopen(lcfname,"w");
      fileOpen=1;
    }
  else return;

  //output header information

  //blending
  strcpy(data,"#fs: ");
  for(int i=0;i<Paramfile->numobservatories;i++)
    {
      sprintf(tmp,"%g ",Event->fs[i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Source magnitudes
  strcpy(data,"#Sourcemag: ");
  int sn = Event->source;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      sprintf(tmp,"%g ",Sources->mags[sn][i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Source data
  strcpy(data,"#Sourcedata: ");
  sprintf(tmp,"%d ",sn); strcat(data,tmp);
  for(int i=0;i<sOutputCols;i++)
    {
      sprintf(tmp,"%g ", Sources->data[sn][sOutputColumns[i]]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Lens magnitudes
  strcpy(data,"#Lensmag: ");
  int ln = Event->lens;
  for(int i=0;i<Paramfile->Nfilters;i++)
    {
      sprintf(tmp,"%g ",Lenses->mags[ln][i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Lens data
  strcpy(data,"#Lensdata: ");
  sprintf(tmp,"%d ",ln); strcat(data,tmp);
  for(int i=0;i<lOutputCols;i++)
    {
      sprintf(tmp,"%g ", Lenses->data[ln][lOutputColumns[i]]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Planet data
  strcpy(data,"#Planet: ");
  for(int i=0;i<NPLANETINPUT+NPLANETDERIV;i++)
    {
      sprintf(tmp,"%g ",Event->params[i]);
      strcat(data,tmp);
    }
  fprintf(lcfile_ptr,"%s\n",data);

  //Microlensing data
  strcpy(data,"#Event: ");
  sprintf(tmp,"%g ",Event->u0); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->alpha); strcat(data,tmp);
  sprintf(tmp,"%.12g ",Event->t0); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->tE_r); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->rs); strcat(data,tmp);
  sprintf(tmp,"%d ",Event->flag_needFS); strcat(data,tmp);
  sprintf(tmp,"%g ",Event->flatchi2); strcat(data,tmp);
  sprintf(tmp,"%g ", 
	  (Event->flag_needFS?Event->FSPL.chisq:Event->PSPL.chisq)); 
  strcat(data,tmp);
  fprintf(lcfile_ptr,"%s\n",data);

  //output the lightcurve

  double mag, magerr, tmag, tmagerr, fs;
  double fs0=Event->fs[0];

  if(lcfile_ptr!=NULL && fileOpen==1)
    {
      for(i=0;i<Event->nepochs;i++)
	{
	  t=Event->epoch[i];
	  obsidx=Event->obsidx[i];
	  fs=Event->fs[obsidx];

	  tmag = -2.5*log10(fs0*(Event->Atrue[i]-1.0+fs)/fs + 1-fs0);
	  mag = - 2.5*log10(fs0*(Event->Aobs[i]-1.0+fs)/fs + 1-fs0);
	  magerr = 1.085736205 * Event->Aerr[i]/Event->Aobs[i] * fs0/fs;
	  tmagerr = 1.085736205 * Event->Atrueerr[i]/Event->Atrue[i] * fs0/fs;
	  
	  fprintf(lcfile_ptr, "%.11g %.8g %g %.8g %g %d %d\n",Event->epoch[i], mag, magerr, tmag, tmagerr, obsidx, (Event->nosat[i]?0:1)); 


	}

      //calculate the theoretical lightcurve

      for(double t=Event->t0-10.0*Event->tE_r; t<=Event->t0+10.0*Event->tE_r;t+=0.01)
      	{
      	  double thAmp = (1.0-Event->fs[0]) + Event->fs[0]*wittFSMagnification(sqrt(sqr(Event->u0)+sqr((t-Event->t0)/Event->tE_r)),Event->rs);
      	  fprintf(lcfile_ptr, "%.11g 1 1 %.8g 1 99 0\n",t, -2.5*log10(thAmp));

      	}
      
      fclose(lcfile_ptr);
    }
}

void outputImages(struct event *Event, struct obsfilekeywords World[], struct slcat* Sources, struct filekeywords* Paramfile)
{
  char tmp1[1000];
  string basefname;
  string extension;
  string oname;
  string imtype;
  int filter;
  double mag;
  double peaktime;
  double background;

  if(!Paramfile->outputImages
     || !(Paramfile->outputOnErr || Paramfile->outputOnDet 
	  || Paramfile->outputOnAll)) return;
  else
    {
      if((Event->lcerror+Event->deterror))
	{
	  if(!Paramfile->outputOnErr) return;
	}
      else if((Event->detected))
	{
	  if(!Paramfile->outputOnDet) return;
	}
      else if((!Event->detected))
	{
	  if(!Paramfile->outputOnAll) return;
	}
    }

  if(Event->detected) extension=string(".det");
  else if(Event->lcerror||Event->deterror) extension=string(".err");
  else extension=string(".all");

  sprintf(tmp1,"%s%s_%d_%d",Paramfile->outputdir, Paramfile->run_name,
	  Event->instance,  Event->id);
  basefname=tmp1;

  //output test images
  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      sprintf(tmp1,".%d_",obsidx);

      filter = World[obsidx].filter;

      //first the baseline image
      imtype=string("base");
      oname = basefname + tmp1 + imtype + extension + string(".fits");

      mag = Sources->mags[Event->source][filter];

      //calculate the background at random time - may as well be first point
      background = Event->backmag[Event->nepochsvec[obsidx]];
      World[obsidx].im.set_background(background);
      World[obsidx].im.addbg();

      //add the baseline source
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], mag);
      World[obsidx].im.reset_detector();
      World[obsidx].im.expose(World[obsidx].exptime[0], 
			      World[obsidx].nstack[0]);
      World[obsidx].im.write_fits(oname,true);
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], mag);
      World[obsidx].im.subbg();

      //last the peak image
      imtype=string("peak");
      oname = basefname + tmp1 + imtype + extension + string(".fits");

      if(Event->Amax < 0) 
	{
	  //if the peak wasn't recorded
	  Event->Amax = (sqr(Event->u0)+2.0)/sqrt(sqr(Event->u0)*(sqr(Event->u0)+4.0));
	  peaktime = Event->t0;
	}
      else
	{
	  peaktime = Event->epoch[Event->peakpoint];
	}

      //background
      World[obsidx].im.set_background(Event->backmag[Event->peakpoint]);
      World[obsidx].im.addbg();

      //add the peak source    
      World[obsidx].im.addstar(Event->xsub[obsidx], Event->ysub[obsidx], 
				   mag-2.5*log10(Event->Amax));
      World[obsidx].im.reset_detector();
      World[obsidx].im.expose(World[obsidx].exptime[0], 
			      World[obsidx].nstack[0]);
      World[obsidx].im.write_fits(oname,true);
      World[obsidx].im.substar(Event->xsub[obsidx], Event->ysub[obsidx], 
			       mag-2.5*log10(Event->Amax));
      World[obsidx].im.subbg();
      //cout << "Amax = " << Event->Amax << endl;
    }
}
