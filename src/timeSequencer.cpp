#include<string>
#include<iostream>
#include <iomanip>

#include "timeSequencer.h"
#include "zodiacalLight.h"
#include "ephem.h"

#define DEBUGVAR 0

void computeAltitude(double JD, double obslong, double obslat, double eventRA, double eventDEC, int obsidx, double *Alt, double *Az);

int inboundary(double x, double y, double X[], double Y[], int ndata);
void whichFields(struct obsfilekeywords World[], struct event *Event, int nobs);
void determineEpochs(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile);
//void computeParallax(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile);
void setupMemory(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile);
void clearVectors(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile);

void timeSequencer(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile, struct slcat* Lenses, struct slcat* Sources)
{
  //Clear old vectors
  clearVectors(World,Event,Paramfile);

  /* Convert Galactic l,b of event to RA and DEC */
  if(DEBUGVAR) cout << "eq2gal" << endl;
  eq2gal(Event->l*TO_RAD,Event->b*TO_RAD,'g',&Event->ra,&Event->dec); 
  if(Event->ra<0) Event->ra+=2*PI; 
  
  /* printf("ra: %f dec: %f\n",Event->ra,Event->dec); */  
  

  /* Find which fields the event is seen by for each observatory */
  if(DEBUGVAR) cout << "whichFields" << endl;
  whichFields(World,Event, Paramfile->numobservatories);     

  /* Show which fields the event is seen by for each observatory */
  /* showIsseenby(World,Event, Paramfile->numobservatories); */      

  /* Find observatory epochs which satisfy weather, moon and altitude conditions. Load Event epoch array */
  if(DEBUGVAR) cout << "determineEpochs" << endl;
  determineEpochs(World,Event,Paramfile);

  /* Calculate parallax shifts */
  //if(DEBUGVAR) cout << "computeParallax" << endl;
  //computeParallax(World,Event,Paramfile);

  //Set up memory
  if(DEBUGVAR) cout << "setupMemory" << endl;
  setupMemory(World,Event,Paramfile);
  
}

void determineEpochs(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile)
{
  /* Altitude limit to observe THIS SHOULD BE IN PARAMETERFILE  */
  double EVENT_ALT_LIM = 20*(TO_RAD);    
  /* Extinction coef. in V. THIS SHOULD BE IN PARAMETERFILE */
  double C_EXT = 0.3;  //set in observatory loop
  /* Sky brightness in V. THIS SHOULD BE IN PARAMETERFILE */
  double VSKY = 21.7;                    

  int obsidx;   /* loop through all observatories */
  int ind;      /* loop through observatory arrays */
  int idx=0;    /* This increments the final epoch array */

  int jnd;      /*This loops through fields*/
  int thisseen; /*Is this event seen by the observatory*/

  double Alt,Az,DeltaV, D ,ObjMoonDist,K;

  bool observing;

  /* Loop through all observatories */
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      thisseen=0;

      C_EXT = World[obsidx].extcoeff;
      VSKY =  World[obsidx].skybackground; //20.0 - 2.5*log10(World[obsidx].constbackground);

      /* Loop through all fields for observatory obsidx */   
      for(jnd=0;jnd<World[obsidx].nfields;jnd++)
	{
	  

	  if(Event->isseenby[obsidx][jnd]){
	    /* Observed in field jnd in observatory obsidx? */
	    /* printf("Is seen by obs %d in field %d",obsidx,jnd); */
	    //leaving thisseen=1?
	    thisseen=1;
	  
	    /* Loop through all observatory epochs */
	    for(ind=0;ind<World[obsidx].nepochs;ind++)
	      {
		/* Is this epoch for the current field?
		   Did weather permit? */
		   
		if(World[obsidx].field[ind]==jnd &&
		   World[obsidx].weather[ind])
		  {
		    /* Is this a ground based observatory? */
		    if(World[obsidx].space==0)
		      {
			//Ground based is no longer working - needs updating
			/* Compute event alititude */			  

			computeAltitude(World[obsidx].jd[ind],
					World[obsidx].longitude,
					World[obsidx].latitude, 
					Event->ra, Event->dec, obsidx,
					&Alt, &Az); 
	   
			moon_sky_brightness(World[obsidx].jd[ind],
					    Event->ra,Event->dec,
					    World[obsidx].longitude*TO_RAD,
					    World[obsidx].latitude*TO_RAD,
					    World[obsidx].altitude,
					    C_EXT , VSKY,  &DeltaV, &D,
					    &ObjMoonDist, &K);

			if(Paramfile->verbosity>=3)
			  {
			    cout << "ra, dec " << Event->ra*TO_DEG << " " << Event->dec*TO_DEG << endl;
			    cout << " " << World[obsidx].epoch[ind] << " " << setw(14) << fixed << World[obsidx].jd[ind];
			    cout.setf(ios_base::fmtflags(0), ios::floatfield);
			    cout << " " << Alt*TO_DEG << " " << Az*TO_DEG << " " << DeltaV << endl;		    
			  }

			/*Test event altitude ALSO SHOULD TEST MOON ANGLE*/
			if(Alt>EVENT_ALT_LIM) observing=true;
			else observing=false;
		      }
		    else observing=true;
	   
		    if(observing)
		      {

			double skyextcorr = 1.0;
			if(!World[obsidx].space)
			  {
			    //Reduce exposure time and increase background to 
			    //simulate the effect of sky extinction. Increase 
			    //the sky background by the same factor to match 
			    //the reduced exposure time

			    skyextcorr = pow(10,-0.4*C_EXT/cos(piO2-Alt));
			  }

			
			/* Epoch of observation */
			Event->epoch.push_back(World[obsidx].epoch[ind]);
			Event->jdepoch.push_back(ind);
			/* Observatory index */
			Event->obsidx.push_back(obsidx);
			/* Exposure time */
			Event->texp.push_back(World[obsidx].exptime[ind] 
					      * skyextcorr);
			Event->nstack.push_back(World[obsidx].nstack[ind]);

			//backgrounds
			Event->backmag.push_back(20.0
						 -2.5*log10(World[obsidx].constbackground 			  
							    + World[obsidx].zodiflux[ind]
							    + (World[obsidx].space==0?pow(10,-0.4*(World[obsidx].skybackground+DeltaV-20))/skyextcorr:0))
						 );

			//parallax shifts
			
			
			//Ground specific quantities
			if(World[obsidx].space==0)
			  {
			    /* Altitude of observation */
			    Event->alt.push_back(Alt);
			    /* Object-Moon distance, radians */
			    Event->moonObjDist.push_back(ObjMoonDist);
			    /* The change in the V-band sky brightness
			       caused by moonlight. */  
			    Event->deltaVmoon.push_back(DeltaV);
			  }  /* end test altitude limit */

			idx++;
		      }
	    
		  } /* end if in field and good weather */
		
	      }  /* end loop through all observatory epochs */

	    
	    }  /* end if isseenby */

	}  /* end loop through all fields for observatory obsidx */

      /*Cumulative ndata array, 1-indexed*/
      Event->nepochsvec[obsidx+1] = idx;

      if(Paramfile->verbosity>2)
	{
	  cout << "nepochsvec[" << obsidx << "] = " << Event->nepochsvec[obsidx] << " " << Event->nepochsvec[obsidx+1] << endl;
	}

      //World[obsidx].Aseen+=thisseen*Event->eventArea;
      //World[obsidx].Nseen+=thisseen;
      //World[obsidx].Aoccured+=Event->eventArea;
      //World[obsidx].Noccured++;

    }  /* end loop through all observatories */

  Event->nepochs = idx;
  Event->nepochsvec[0] = 0;

  /*Just so that we can pass this number into the fitting algorithm*/
  Event->numobservatories = Paramfile->numobservatories;
}

void whichFields(struct obsfilekeywords World[], struct event *Event, int nobs)
{
  int ind,jnd;
  double Fr[5];
  double Fd[5];

  int obsidx;

  /*  LOOP THROUGH ALL OBSERVATORIES */
  for(obsidx = 0;obsidx<nobs;obsidx++)
    {

      /*  LOOP THROUGH ALL FIELDS  */
      for(ind = 0;ind<World[obsidx].nfields;ind++)
	{
	  /* LOOP THROUGH EACH VERTEX */
	  for(jnd=0;jnd<4;jnd++)
	    {  
	      Fr[jnd] = World[obsidx].fieldVerticies[0][jnd][ind];
	      Fd[jnd] = World[obsidx].fieldVerticies[1][jnd][ind];
	    }
	  
	  /* add if for USE_FIELDS*/
	  Event->isseenby[obsidx][ind] = inboundary(Event->l*TO_RAD,
						    Event->b*TO_RAD, Fr, 
						    Fd, 4);
	}  /* end loop through all fields */
    }  /* end loop through all observatories */
}


int inboundary(double x, double y, double X[], double Y[], int ndata)
{
  /* Find if point is in a given boundary */
  int ind;
  int As=0;
  double A;
  
  X[ndata] = X[0];  
  Y[ndata] = Y[0]; 
 
  for(ind=0;ind<ndata;ind++)
    {
    
      A = 0.5*(x*(Y[ind] - Y[ind+1]) - y*(X[ind] - X[ind+1]) + X[ind]*Y[ind+1] 
	       - Y[ind]*X[ind+1]);
      /* printf("%f %f %f %f %f %d\n",X[ind], X[ind+1], Y[ind], Y[ind+1],A,dsgn(A)); */   
      As+=dsgn(A);
    }
  /*  printf("ndata = %d\n",ndata); */
  
  /* sum of signs of areas should equal number of boundary points for point in
     boundary */
  if(abs(As)==ndata) return(1);
  else return(0);   //return(1);  //FOR TESTING ONLY!!!!!!!!
  
}


void computeAltitude(double JD, double obslong, double obslat, double eventRA, double eventDEC, int obsidx, double *Alt, double *Az)
{
  double LST,alt,az;

  lst(JD,obslong*TO_RAD,&LST);
  eq2horiz(eventRA,eventDEC,JD,obslat*TO_RAD,obslong*TO_RAD,'h',&az,&alt);
  *Alt = alt;
  *Az = az;
}

//Resize the vectors that will hold the lightcurves
void setupMemory(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile)
{

  Event->Aobs.resize(Event->nepochs);
  Event->Aerr.resize(Event->nepochs);
  Event->Atrue.resize(Event->nepochs);
  Event->Atrueerr.resize(Event->nepochs);
  Event->Afit.resize(Event->nepochs);
  Event->nosat.resize(Event->nepochs);
  Event->xs.resize(Event->nepochs);
  Event->ys.resize(Event->nepochs);
  Event->xl1.resize(Event->nepochs);
  Event->yl1.resize(Event->nepochs);
  Event->xl2.resize(Event->nepochs);
  Event->yl2.resize(Event->nepochs);

}

//clear the vectors ready to be filled again
void clearVectors(struct obsfilekeywords World[], struct event *Event, struct filekeywords *Paramfile)
{
  Event->epoch.clear();
  Event->jdepoch.clear();
  Event->alt.clear();
  Event->obsidx.clear();
  Event->nstack.clear();
  Event->texp.clear();
  Event->moonObjDist.clear();
  Event->deltaVmoon.clear();
  Event->nosat.clear();      /*Is point unsaturated? */
  Event->backmag.clear();

}
