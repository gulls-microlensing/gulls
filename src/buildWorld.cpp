/*! \file 
\brief Functions to create the world of observatories

This file contains functions that populate the world with the observatories
listed in the file observatory.list and computes the possible observation times
of the event for each observatory.
 */
#include "buildWorld.h"
#include "strfns.h"
#include "readObservatoryfile.h"
#include "constdefs.h"
#include "ephem.h"
#include "split.h"

#include<gsl/gsl_sf_erf.h>

#include<string>
#include<vector>
#include<fstream>
#include<iostream>

using namespace std;

#define DEBUGVAR 1

void loadFieldCentres(struct obsfilekeywords World[], struct filekeywords *Paramfile, int obsidx);
void buildWorldErr(int err, int line, const char msg[]);
void computeFieldVerticies(struct obsfilekeywords World[], int obsidx, ofstream& logfile_ptr);
void showVerticies(struct obsfilekeywords World[], int obsidx, ofstream& logfile_ptr);
void surveyArea(struct obsfilekeywords World[],int obsidx,ofstream& logfile_ptr);
int inboundary(double x, double y, double X[], double Y[], int ndata);
void loadObsSequence(struct obsfilekeywords* World, struct filekeywords *Paramfile);
void applyObsSequence(struct obsfilekeywords World[], struct filekeywords *Paramfile);
void setupImage(struct obsfilekeywords World[], struct filekeywords *Paramfile, long* idum);
double max(double x, double y);
void setEpochProperties(struct obsfilekeywords World[], struct filekeywords *Paramfile);
void setupOrbit(struct obsfilekeywords World[], struct filekeywords *Paramfile);

void buildWorld(struct filekeywords *Paramfile, struct obsfilekeywords World[], long *idum, ofstream& logfile_ptr)
{
  
  void findSunRiseSetTimes(struct filekeywords *Paramfile, struct obsfilekeywords World[], int numobservatories);
  void applyWeather(struct obsfilekeywords World[], struct filekeywords *Paramfile, long *idum);
  void buildFields(struct obsfilekeywords World[], struct filekeywords *Paramfile, ofstream& logfile_ptr);
  void collectingArea(struct obsfilekeywords World[], int numobservatories);
  char listfile[1000];
  FILE* obslistfile_ptr; 
  char str[100];
  char tmp[100];
  char obsfile[1000];
  int obsctr = 0;
  int allspace=1;
  int i;

  strcpy(listfile,Paramfile->obslist);

  /* Open up listfile containing observatories */
  if(strlen(listfile))       
    { 	
      obslistfile_ptr = fopen(listfile, "r");
      if (obslistfile_ptr == NULL)
	{
	  sprintf(str,"Unable to open file: %s",Paramfile->obslist);
	  fmtline(str,WIDTH,"FAILED");
	  exit(1);
	}
      sprintf(str,"Observatory list file: %s",Paramfile->obslist);
      fmtline(str,WIDTH,"READY");
    }
  else
    {
      sprintf(str,"Unable to find file: %s",Paramfile->obslist);
      fmtline(str,WIDTH,"FAILED");
      exit(1);
    }

  /* For each observatory listed in the observatory list file, read in
     keywords */
  obsctr=0; 
  while (fgets(tmp,sizeof(tmp),obslistfile_ptr) != NULL)
    {
      /* Ignore comment lines and empty lines */
      if(tmp[0]=='#' || strlen(tmp)<=1) continue;

      sprintf(obsfile,"%s%s",Paramfile->obsdir,tmp);
      obsfile[strlen(obsfile)-1]='\0';

      if(DEBUGVAR)  {printf("readObservatoryfile\n"); fflush(stdout);}
      readObservatoryfile(obsfile,World,obsctr);
      if(DEBUGVAR)  {printf("Observatoryfile read\n"); fflush(stdout);}

      sprintf(str,"%s",obsfile); 
      fmtline(str,2*WIDTH,"PARSED"); 
      World[obsctr].Aseen=0;
      World[obsctr].Aoccured=0;
      World[obsctr].Nseen=0;
      World[obsctr].Noccured=0;

      obsctr++;
    }

 
  Paramfile->numobservatories = obsctr;
  
  /*check if all observatories are space observatories*/
  for(i=0;i<Paramfile->numobservatories;i++)
    {
      if(!World[i].space) allspace=0;
    }

  /* For each observatory, compute collecting area */
  if(DEBUGVAR)  {printf("collectingArea\n"); fflush(stdout);}
  collectingArea(World, Paramfile->numobservatories);

  /* For each observatory compute sun rise/set times */
  if(DEBUGVAR)  {printf("findSunRiseSetTimes\n"); fflush(stdout);}
  if(!allspace)
    {
      findSunRiseSetTimes(Paramfile, World,Paramfile->numobservatories);
    }

  /* For each observatory, load its observing sequence */
  if(DEBUGVAR) {printf("loadObsSequence\n"); fflush(stdout);}
  loadObsSequence(World, Paramfile);

  if(DEBUGVAR)
    {
      int onx, inx;
      logfile_ptr << "Observing sequence\n";
      for(onx=0;onx<Paramfile->numobservatories;onx++)
	{
	  logfile_ptr << "\nObservatory " << onx << ":\n";
	  for(inx=0;inx<World[onx].sequence_length;inx++)
	    {
	      logfile_ptr << World[onx].sequence[inx].field << " " 
			   << World[onx].sequence[inx].nstack << " "
			   << World[onx].sequence[inx].texp << endl;
	    }
	}
    }

  /* For each observatory compute instrument cadence and populate time
     vectors */
  if(DEBUGVAR)  {printf("applyObsSequence\n"); fflush(stdout);}
  applyObsSequence(World,Paramfile);

  /* For each observatory apply weather profiles*/
  if(DEBUGVAR)  {printf("applyWeather\n"); fflush(stdout);}
  applyWeather(World,Paramfile, idum);
 
  /*For each observatory, load field centres, compute field vertices */
  if(DEBUGVAR)  {printf("buildFields\n"); fflush(stdout);}
  buildFields(World,Paramfile,logfile_ptr);

  if(DEBUGVAR)  {printf("setEpochProperties\n"); fflush(stdout);}
  setEpochProperties(World, Paramfile);

  if(DEBUGVAR) {printf("setupImage\n"); fflush(stdout);}
  setupImage(World,Paramfile,idum);

  if(DEBUGVAR) {printf("setupOrbit\n"); fflush(stdout);}
  setupOrbit(World,Paramfile);

  sprintf(str,"World created: %d observatories",Paramfile->numobservatories);
  fmtline(str,WIDTH,"READY");

}


/*! Compute collecting area for each telescope */
void collectingArea(struct obsfilekeywords World[], int numobservatories)
{
  int obsidx;
  for(obsidx = 0;obsidx<numobservatories;obsidx++)
    {
      World[obsidx].collectingarea = PI*(pow(0.5*World[obsidx].primary,2) 
					 - pow(0.5*World[obsidx].blockage,2));

    }

}

/*! Compute field verticies using the set of field centres for each observatory */
void buildFields(struct obsfilekeywords World[], struct filekeywords *Paramfile, ofstream& logfile_ptr)
{ 
  int obsidx;
  
  for(obsidx = 0;obsidx<Paramfile->numobservatories;obsidx++)
    {

      loadFieldCentres(World,Paramfile,obsidx);
      computeFieldVerticies(World,obsidx,logfile_ptr);
      
      showVerticies(World, obsidx, logfile_ptr); 
    
      surveyArea(World,obsidx,logfile_ptr);

    }
}


/*! Find the sunrise and sunset times for each observatory for the simulation */
void findSunRiseSetTimes(struct filekeywords *Paramfile, struct obsfilekeywords World[], int numobservatories)
{
  int MAX_NUM_ALT = Paramfile->NUM_SIM_DAYS*24*60/ALT_STEP;

  double Ra,Dec,R,SL,EquationTime,Alt,Az;
  vector<double> ALT(MAX_NUM_ALT); 
  vector<double> TVEC(MAX_NUM_ALT); 
  double tincr,d0,d1;
  int ind,signumsum,idx0,idx1,fdx;
  int start=0;
  double ALTLIM = -6 * (TO_RAD); /* ALTITUDE LIMIT DEFINED HERE  */
  int obsidx;

  double sim_t0 = Paramfile->simulation_zerotime;

  /* Compute altitude for whole simulation time. */
  /* printf("tinc\n"); */

  /* time step for whole simulation time vector */
  tincr =  (double)ALT_STEP / 1440.0; 

  /* printf("for each obs\n"); */
  for(obsidx=0;obsidx< numobservatories;obsidx++)
    {

      World[obsidx].sunriseset = vector<vector<double> >(Paramfile->NUM_SIM_DAYS+4,vector<double>(2));

      /* No need to compute sunrise for space telescope */
      if(World[obsidx].space==1) continue; 
      
      for(ind=0;ind<MAX_NUM_ALT;ind++)
	{
	  TVEC[ind] = ind*tincr;
	  
	  suncoo(TVEC[ind]+sim_t0, 'j', &Ra, &Dec, &R,
		 &SL, &EquationTime);
	  eq2horiz(Ra, Dec, TVEC[ind]+sim_t0, 
		   World[obsidx].latitude*(TO_RAD), 
		   World[obsidx].longitude*(TO_RAD), 'h', &Az,&Alt);
	  if(Paramfile->verbosity>3)
	    cout << TVEC[ind] << " " << TVEC[ind] + sim_t0 << " " << Az << " " << Alt << endl;
	  ALT[ind] =  Alt;
	}
 
      /* print2Vec(TVEC,ALT,MAX_NUM_ALT,"altitude.dat");  */
      if(Paramfile->verbosity>2)
	cout << "sunriseset size " << World[obsidx].sunriseset.size() << endl;
      if(Paramfile->verbosity>2)
	cout << "sunriseset size2 " << World[obsidx].sunriseset[0].size() << endl;

      idx1=0; idx0=0; start=0;

      for(ind=0;ind<MAX_NUM_ALT-1;ind++)
	{
	  d0 =  ALT[ind] - ALTLIM;
	  d1 = ALT[ind+1] - ALTLIM;

	  signumsum = dsgn(d0) + dsgn(d1);

	  /* Did consecutive altitudes cross altitude limit ALTLIM? */
	  if(signumsum ==0)
	    { 

	      /* NB: WE ARE ALWAYS TAKING THE LATEST SET TIME AND EARLIEST
		 RISE TIME */
	      if( ALT[ind+1] > ALTLIM)
		{
		  if(start) World[obsidx].sunriseset[idx0][0] = TVEC[ind];
		  /*if(start)  printf("Rise: %f %f %f %f\n",
		    ALT[ind]*(TO_DEG),ALT[ind+1]*(TO_DEG),TVEC[ind],
		    TVEC[ind+1]); */
		  if(start)  idx0++;
		}
    
	      if(ALT[ind+1] < ALTLIM)
		{
		  /* This is to ensure that we start from the first sunset */
		  start = 1; 
		  /* printf("Set: %f %f %f %f\n",
		     ALT[ind]*(TO_DEG),ALT[ind+1]*(TO_DEG),
		     TVEC[ind],TVEC[ind+1]); */
		  World[obsidx].sunriseset[idx1][1] = TVEC[ind+1];
		  idx1++;
		}
	    }
	}

      /* This is to ensure that we have both rise and set times */
      if(idx0<idx1) fdx = idx0; 
      else fdx = idx1;

      /* printf("final idx\n"); */
      World[obsidx].finalidxnights = fdx;

      /* for(ind=0;ind<World[obsidx].finalidxnights;ind++) 
	 printf("%d %f %f\n",obsidx, World[obsidx].sunriseset[ind][0],
	 World[obsidx].sunriseset[ind][1]);    */
      /* printf("done\n"); */

    }
}

/*! Compute the set of all possible observations for each observatory */
void applyObsSequence(struct obsfilekeywords World[], struct filekeywords *Paramfile)
{ 
  int obsidx;
  int ind = 0;
  int jdx,ldx;
  double tx,tend;
  int lastfield=-1;
  int field,nstack;
  double texp;
  double tread;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //time taken to read CCD array
      tread = World[obsidx].readohead/SECINDAY;
      //cout << "obs " << obsidx << " " << tread << endl;

      if(World[obsidx].space==0)
	{
	  ind=0;
	  
	  if(DEBUGVAR) 
	    {
	      printf("Observatory %d Sequence length: %d\n",obsidx,
		     World[obsidx].sequence_length);
	      fflush(stdout);
	    }
	  
	  //Loop through sunsetrise times
	  for(jdx=0;jdx<World[obsidx].finalidxnights;jdx++)
	    {
	      
	      //Set time to (jdx)th sunset
	      tx = World[obsidx].sunriseset[jdx][1];
	      
	      while (tx<World[obsidx].sunriseset[jdx][0])
		{
		  for(ldx=0;ldx<World[obsidx].sequence_length;ldx++)
		    {
		      field = World[obsidx].sequence[ldx].field;
		      nstack = World[obsidx].sequence[ldx].nstack;
		      texp = World[obsidx].sequence[ldx].texp/SECINDAY;

		      //Slewing and readout etc are handled by the observation
		      //sequence files, unless it is the readout between
		      //images of stacks, which gets added on here
		      //To denote readout or slewing use Nstack=-1 
		      //Texp=-tslew/tread

		      //Stop if the sun has risen and park the scope
		      if(tx>=World[obsidx].sunriseset[jdx][0]) 
			{
			  lastfield = -1;
			  break;
			}

		      //Were we actually observing (i.e. no if using another
                      //instrument)
		      if(nstack>0 && texp>0)
			{
			  World[obsidx].epoch.push_back(tx + 0.5*(texp*nstack + tread*(nstack-1)));
			  World[obsidx].field.push_back(field);
			  World[obsidx].nstack.push_back(nstack);
			  World[obsidx].exptime.push_back(texp*SECINDAY);
			  ind++;
			  if (ind>=MAX_NUM_EPOCH)
			    buildWorldErr(1,__LINE__,"");
			}

		      tx += abs(nstack)*abs(texp) + (abs(nstack)-1)*tread;
		      //the final readout time is handled by the sequence
		      
		    } /* end loop through observation sequence */

		} /* end while sun is down */
	      
	    } /* end loop through sunrise times */
	  
	} /* end if ground based */

      if(World[obsidx].space==1)
	{
	  ind=0;
	  tx = 0; //zerotime is added only when needed
	  tend = Paramfile->NUM_SIM_DAYS;

	  if(DEBUGVAR) 
	    {
	      printf("Observatory %d Sequence length: %d\n",obsidx,
		     World[obsidx].sequence_length);
	      fflush(stdout);
	    }
      
	  while (tx<tend)
	    {
	      for(ldx=0;ldx<World[obsidx].sequence_length;ldx++)
		{
		  field=World[obsidx].sequence[ldx].field;
		  nstack=World[obsidx].sequence[ldx].nstack;
		  texp=World[obsidx].sequence[ldx].texp/SECINDAY;

		  //Slew if moving to a new field we actually observe
		  if(field!=lastfield && field>=0 
		     && field<World[obsidx].nfields)
		    {
		      lastfield=field;
		    }

		  //cout << field << " " << texp << " " << nstack << " " << tx << endl;

		  //Were we actually observing (i.e. no if using another
		  //instrument)
		  if(nstack>0 && texp>0)
		    {
		      World[obsidx].epoch.push_back(tx + 0.5*(texp*nstack + tread*(nstack-1)));
		      World[obsidx].field.push_back(field);
		      World[obsidx].nstack.push_back(nstack);
		      World[obsidx].exptime.push_back(texp*SECINDAY);
		      ind++;
		      if (ind>=MAX_NUM_EPOCH)
			buildWorldErr(1,__LINE__,"");
		    }

		  tx += abs(nstack)*abs(texp) + (abs(nstack)-1)*tread;
		  //the final readout time is handled by the sequence
		  
		} /* end loop through observation sequence */
	      
	    } /* end while we are in simulation time */

	  if (ind>=MAX_NUM_EPOCH)
	    {
	      printf("ind=%d; MNE=%d\n",ind,MAX_NUM_EPOCH);
	      buildWorldErr(1,__LINE__,"");
	      fflush(stdout);
	    }
	}

      World[obsidx].nepochs=ind;

    } /* end loop over observatories */

}

//Calculate various properties that apply universally to each epoch
void setEpochProperties(struct obsfilekeywords World[], struct filekeywords *Paramfile)
{
  int idx, obsidx;
  double jd;
  int field;
  double lmlsun;

  char throughput[1000];

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      //setup zodi
      strcpy(throughput,Paramfile->obsdir);
      strcat(throughput,World[obsidx].throughput);
      World[obsidx].zodi.set_bandpass(string(throughput));

      for(idx=0;idx<World[obsidx].nepochs;idx++)
	{
	  field = World[obsidx].field[idx];

	  //Julian date
	  World[obsidx].jd.push_back(jd = World[obsidx].epoch[idx]+Paramfile->simulation_zerotime);

	  //Solar position
	  World[obsidx].lambdasun.push_back(solarlambda(jd));

	  World[obsidx].lambda.push_back(0);
	  World[obsidx].beta.push_back(0);

	  //field position in ecliptic coordinates
	  gal2eclip(World[obsidx].fieldCentres[0][field], World[obsidx].fieldCentres[1][field], &World[obsidx].lambda[idx], &World[obsidx].beta[idx]);
 
	  
	  lmlsun=World[obsidx].lambda[idx]-World[obsidx].lambdasun[idx];	  
	  while(lmlsun<0) lmlsun+=360;

	  //Zodiacal light level - in mag20 per sq arcsec
	  World[obsidx].zodiflux.push_back(pow(10,-0.4*(World[obsidx].zodi.get_mag(lmlsun,World[obsidx].beta[idx])-20)));
	  //World[obsidx].zodiflux.push_back(0);
	  //cout << "zodiflux = " << World[obsidx].zodiflux.back() << " " << World[obsidx].zodi.get_mag(lmlsun,World[obsidx].beta[idx]) << " " << lmlsun << " " << World[obsidx].beta[idx]<< endl;
	}
    }
}


/*! Apply the weather profile for each observatory to the set of possible observations 
  computed in computeCadence */
void applyWeather(struct obsfilekeywords World[], struct filekeywords *Paramfile, long *idum)
{
  void loadWeatherProfile(char profileName[], vector<double>* WeatherProfile, struct filekeywords *Paramfile);
  int getWeather(double PrClearnight, long *idum);
  void dayofyear(int day, int month, int year, int *doy);
  
  int obsidx;
  int ind;
  int weatherval=0;
  int currentday;  //number of quarter days since beginning of simulation
  int oldday=-50000000;

  /* printf("for each obs\n"); */
  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      /* printf("for each epoch; nepochs =  %d\n",World[obsidx].nepochs); */
      if(DEBUGVAR) {printf("loadWeatherProfile %s\n",
			  World[obsidx].weatherProfile); fflush(stdout);}
      loadWeatherProfile(World[obsidx].weatherProfile, 
			 &World[obsidx].weatherSequence, Paramfile);
      if(DEBUGVAR) {printf("profileLoaded\n"); fflush(stdout);}
      /* printf("for each epoch; nepochs =  %d\n",World[obsidx].nepochs); */

      /* loop through World[obsidx].epoch and set World[obsidx].weather */
      for(ind=0;ind<World[obsidx].nepochs;ind++)
	{ 

	  currentday = int(floor(World[obsidx].epoch[ind]*4));

	  if(currentday!=oldday)
	    {
	      weatherval 
		= getWeather(World[obsidx].weatherSequence[currentday], idum);
	      oldday=currentday;
	    }

	  World[obsidx].weather.push_back(weatherval);
	  
	}

    }
} 

/*! Roll dice to see whether the weather today permits observations */
int getWeather(double PrClearnight, long *idum)
{
  if(ran2(idum)<PrClearnight) return(1);
  else return(0);
}

/*! Read in weather profile */
void loadWeatherProfile(char profileName[], vector<double>* WeatherProfile, struct filekeywords *Paramfile)
{
  FILE *infile_ptr;
  char filename[1000];
  strcpy(filename,Paramfile->weatherprofiledir);
  strcat(filename,profileName);

  *WeatherProfile = vector<double>((Paramfile->NUM_SIM_DAYS+1)*4);
  
  /* printf("load weather profile\n"); */
  
  /* printf("%s\n",filename); */

  if(DEBUGVAR) {printf("loadWeather opening: %s\n",filename); fflush(stdout);}
  infile_ptr = fopen(filename,"r");
  int ind = 0;
  double null;
  int nindmax;
  int i, ncleardays=0;
  double tmp;
  
  if (infile_ptr == NULL)  buildWorldErr(2, __LINE__,filename); 
  
  while (fscanf(infile_ptr, "%lf %lf", &null, &tmp) != EOF) 
    {
      (*WeatherProfile)[ind] = tmp;
      ind++;
      if(ind==(Paramfile->NUM_SIM_DAYS+1)*4) break;
    }
  fclose(infile_ptr);

  /* if there are any spare simulation days, fill them with zero */

  nindmax=ind;
  {printf("nindmax=%d; NUM_SIM_DAYS=%d\n",nindmax,Paramfile->NUM_SIM_DAYS); fflush(stdout);}
  while(ind<(Paramfile->NUM_SIM_DAYS+1)*4)
    {
      (*WeatherProfile)[ind]=0;
      ind++;
    }

  if(DEBUGVAR)
    {
      printf("Weather Profile:\n"); fflush(stdout);
      for(i=0;i<ind;i++) 
	{
	  printf("%d", int((*WeatherProfile)[i]));
	  ncleardays+=int((*WeatherProfile)[i]);
	}
      
      printf("\nncleardays=%f\n",ncleardays/4.0); fflush(stdout);
    }


  /* printf("loaded weather profile\n"); */

}

/*! Read in observation sequences of all observatories */
void loadObsSequence(struct obsfilekeywords World[], struct filekeywords *Paramfile)
{
  char filename[1000];
  int obsidx;

  ifstream infile;

  int ind = 0;

  vector<string> data;
  string line;
  int repeating=0;
  int repeatCount=0;
  int linecount=0;

  string repeatCommand = string("BEGIN_REPEAT");
  string endRepeatCommand = string("END_REPEAT");

  obssequence repeatBuffer;
  vector<obssequence> repeatStack;

  for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      strcpy(filename,Paramfile->obsdir);
      if(!Paramfile->identicalSequence)
	strcat(filename,World[obsidx].observationSequence);
      else
	strcat(filename,World[0].observationSequence);

      //check if this is a repeat sequence (i.e. can we reuse the lc)
      /*World[obsidx].same_sequence=-1;
      for(int o=0; o<obsidx; o++)
	{
	  if(World[obsidx].same_sequence<0 &&
	     !strcmp(World[obsidx].observationSequence,
		     World[o].observationSequence))
	    {
	      //World[obsidx].same_sequence=o; //a match!
	      //will still need assign epochs, zodi etc
	    }
	    }*/
      
      if(DEBUGVAR) {printf("loadObsSequence opening: %s\n",filename); fflush(stdout);}
      infile.open(filename);
      
      if (!infile)  
	{
	  buildWorldErr(4, __LINE__,filename); 
	}
      
      ind=0;
      repeating=0;
      repeatCount=0;
      linecount=0;
      repeatStack.clear();
      World[obsidx].mintexp=1e10;
      
      while(!infile.eof())
	{
	  getline(infile,line);
	  
	  if(line.compare(0,1,"#"))
	    {
	      split(line,data);
	      
	      if(data.size()>=3)
		{
		  if(!repeating)
		    {
		      World[obsidx].sequence[ind].field = atoi(data[0].c_str());
		      World[obsidx].sequence[ind].nstack = 
			atoi(data[1].c_str());
		      World[obsidx].sequence[ind].texp = atof(data[2].c_str());

		      if(World[obsidx].sequence[ind].nstack>0&&World[obsidx].sequence[ind].texp>0)
			{
			  if(World[obsidx].sequence[ind].texp
			     <World[obsidx].mintexp)
			    {
			      World[obsidx].mintexp = 
				World[obsidx].sequence[ind].texp;
			    }
			}
		      
		      ind++;
		    }
		  else
		    {
		      repeatBuffer.field = atoi(data[0].c_str());
		      repeatBuffer.nstack = atoi(data[1].c_str());
		      repeatBuffer.texp = atof(data[2].c_str());

		      if(repeatBuffer.nstack>0&&repeatBuffer.texp>0)
			{
			  if(repeatBuffer.texp<World[obsidx].mintexp)
			    World[obsidx].mintexp = repeatBuffer.texp;
			}

		      repeatStack.push_back(repeatBuffer);
		    }
		}
	      else if(data.size()>=2)
		{
		  if(data[0].find(repeatCommand)!=string::npos)
		    {
		      repeatCount=atoi(data[1].c_str());
		      repeating=1;
		    }
		}
	      else if(data.size()>=1)
		{
		  if(data[0].find(endRepeatCommand)!=string::npos)
		    {
		      repeating=0;

		      //output the repeat stack the right amount of times
		      while(repeatCount>0)
			{
			  for(int rep=0; rep<int(repeatStack.size()); rep++)
			    {
			      World[obsidx].sequence[ind] = repeatStack[rep];
			      ind++;

			      if(ind==MAX_SEQUENCE_LENGTH) 
				{
				  printf("Warning: observation sequence length exceeded, simulation may not behave as expected\n");
				  break;
				} //end if ind==MAX
			    } //end for each rep
			  repeatCount--;
			} //end while repeatCount
		      //clear the repeat stack
		      repeatStack.clear();
		    } //end if find end repeat
		} //end if data size >=1

	      if(ind==MAX_SEQUENCE_LENGTH) 
		{
		  printf("Warning: observation sequence length exceeded, simulation may not behave as expected\n");
		  break;
		}
	    }

	  linecount++;

	} //end while read line

      if(Paramfile->verbosity) cout << "mintexp[" << obsidx << "] = " << World[obsidx].mintexp << endl;
      
      infile.close();
      infile.clear();

      //check if an end repeat was forgotten
      if(repeating)
	{
	  repeating=0;
	  
	  //output the repeat stack the right amount of times
	  while(repeatCount>0)
	    {
	      for(int rep=0; rep<int(repeatStack.size()); rep++)
		{
		  World[obsidx].sequence[ind] = repeatStack[rep];
		  ind++;

		  if(ind==MAX_SEQUENCE_LENGTH) 
		    {
		      printf("Warning: observation sequence length exceeded, simulation may not behave as expected\n");
		      break;
		    } //end if ind==MAX
		} //end for each rep
	      repeatCount--;
	    } //end while repeatCount
	  repeatStack.clear();
	} //end if repeating
	
      World[obsidx].sequence_length=ind;
    }

}

/*! Read in field centres*/
void loadFieldCentres(struct obsfilekeywords World[], struct filekeywords *Paramfile, int obsidx)
{
  FILE *infile_ptr;
  char filename[1000];
  char str[200];
  strcpy(filename,Paramfile->obsdir);
  strcat(filename,World[obsidx].fieldCentreFile);
  
  /* printf("%d\n",Paramfile->numobservatories); */
  /*   printf("%s\n",filename); */
  
  /*  printf("%s\n",filename);  */
  /*filename[strlen(filename)-1] = '\0';*/

  infile_ptr = fopen(filename,"r");
  int ind = 0;
 
  if (infile_ptr == NULL)  buildWorldErr(2, __LINE__,filename); 
 
  while (fscanf(infile_ptr, "%lf %lf", &World[obsidx].fieldCentres[0][ind],
		&World[obsidx].fieldCentres[1][ind]) != EOF) 
    {
      if(ind>=MAX_NUM_FIELDS) buildWorldErr(3, __LINE__,filename); 
      ind++;
    }
  fclose(infile_ptr);

  if(ind!=World[obsidx].nfields)
    {
      sprintf(str,"Number of field centres (%d) in",ind);
      fmtline(str,50,"(loadFieldCentres)");
      sprintf(str, "%s\ndoes not equal NFIELDS (%d) in config file for observatory %s (buildWorld:loadFieldCentres)",
	      filename,World[obsidx].nfields, World[obsidx].name);
      fmtline(str,100,"");
      sprintf(str,"Exiting");
      fmtline(str,50,"ERROR");
      exit(1);
    }

}

/*! Compute the field verticies for all fields for given observatory */
void computeFieldVerticies(struct obsfilekeywords World[], int obsidx, ofstream& logfile_ptr)
{
  double B2[4] = {0.5,0.5,-0.5,-0.5};
  double B1[4] = {0.5,-0.5,-0.5,0.5};
  double RA,DEC;
  int ind,jnd;
  
  double eta,xsi,beta,gamma,d,a;

  for (jnd = 0; jnd<World[obsidx].nfields;jnd++)
    {

      RA = World[obsidx].fieldCentres[0][jnd]*TO_RAD;
      DEC = World[obsidx].fieldCentres[1][jnd]*TO_RAD;

      for (ind=0;ind<4;ind++)
	{
	  
	  eta = B1[ind]*TO_RAD*World[obsidx].npixy
	    * World[obsidx].pixelsize/3600.0;
	  xsi = B2[ind]*TO_RAD
	    * World[obsidx].npixx*World[obsidx].pixelsize/3600.0;

	  beta = cos(DEC) - eta * sin(DEC);
	  a = atan2(xsi, beta) + RA;
	  gamma = sqrt((xsi*xsi) +(beta*beta));
	  d = atan2(eta*cos(DEC)+sin(DEC) , gamma);
	  
	  World[obsidx].fieldVerticies[0][ind][jnd] = a;
	  World[obsidx].fieldVerticies[1][ind][jnd] = d;

	}
    }
}


/*! Show the field verticies for all fields for given observatory */
void showVerticies(struct obsfilekeywords World[], int obsidx, ofstream& logfile_ptr)
{
  int ind,jnd;
  double l,b;
  double ra,dec;

  for(jnd=0;jnd<World[obsidx].nfields;jnd++)
    {
      for (ind=0;ind<4;ind++) 
	{
	  ra=World[obsidx].fieldVerticies[0][ind][jnd];
	  dec=World[obsidx].fieldVerticies[1][ind][jnd];
	  eq2gal(ra, dec,'g',&l,&b);
	  /*  fprintf(logfile_ptr,
	      "Observatory: %d Field %d RA: %f DEC: %f L: %f B: %f\n",
	      obsidx,jnd,ra,dec,l,b); */
	  /* Reversed order as fields now specified in galactic coordinates */
	  logfile_ptr << "Observatory: " << obsidx << " Field " << jnd << " RA: " << l << " DEC: " << b << " L: " << ra << " B: " << dec << "\n";
	}   
    }
}

/*! Define error messages for buildWorld function */
void buildWorldErr(int err, int line, const char msg[])
{
  char str[100];
  
  sprintf(str,"Error in buildWorld.c at line: %d", line); 
  fmtline(str,50,"(BuildWorldErr)");

  switch (err)
    {

    case 1:
      sprintf(str,"Array overflow (MAX_NUM_EPOCH)"); 
      fmtline(str,50,"FATAL"); exit(1);
      break;

    case 2:
      sprintf(str,"Unable to open profile file: %s",msg); 
      fmtline(str,50,"FATAL"); exit(1);
      break;

    case 3:
      sprintf(str,"Array overflow (MAX_NUM_FIELDS)"); 
      fmtline(str,50,"FATAL"); exit(1);
      break;

    case 4:
      sprintf(str,"Unable to open obs sequence file: %s",msg);
      fmtline(str,50,"FATAL"); exit(1);

    default:
      break;
      
    }
}

void surveyArea(struct obsfilekeywords World[],int obsidx, ofstream& logfile_ptr)
{
  /*Survey area needs to only measure unique area*/

  double fieldVerticies[3][4*MAX_NUM_FIELDS+1]; /*x,y,angle*/

  int hnd,ind,jnd,knd,i;
  int within;

  int nv;

  double Fr[5];
  double Fd[5];

  double minx,maxx,miny,maxy,midx,midy;

  int swapped;

  double area;

  double l,b;

  double t[3];

  nv=0;
  minx=miny=1.0e30; maxx=maxy=-1.0e30;

  /*So find field vertices not within other fields*/
  for(knd=0;knd<World[obsidx].nfields;knd++)
    {
      for(ind=0;ind<4;ind++)
	{
	  /*printf("Vertices: %f %f\n",
	    World[obsidx].fieldVerticies[0][ind][knd],
	    World[obsidx].fieldVerticies[1][ind][knd]);*/
	  within=0;
	  for(jnd=0;jnd<World[obsidx].nfields;jnd++)
	    {
	      if(jnd!=knd)
		{	        
		  /* LOOP THROUGH EACH VERTEX */
		  for(hnd=0;hnd<4;hnd++)
		    {
		      Fr[hnd] = World[obsidx].fieldVerticies[0][hnd][jnd];
		      Fd[hnd] = World[obsidx].fieldVerticies[1][hnd][jnd];
		    }
		
		  within+=inboundary(World[obsidx].fieldVerticies[0][ind][knd],
				     World[obsidx].fieldVerticies[1][ind][knd],
				     Fr,Fd,4);
		}
	    }

	  if(!within)
	    {
	      fieldVerticies[0][nv]=World[obsidx].fieldVerticies[0][ind][knd];
	      /*if(fieldVerticies[0][nv]>PI) fieldVerticies[0][nv]-=2.0*PI;*/
	      if(fieldVerticies[0][nv]>maxx) maxx=fieldVerticies[0][nv];
	      if(fieldVerticies[0][nv]<minx) minx=fieldVerticies[0][nv];
	      fieldVerticies[1][nv]=World[obsidx].fieldVerticies[1][ind][knd];
	      if(fieldVerticies[1][nv]>maxy) maxy=fieldVerticies[1][nv];
	      if(fieldVerticies[1][nv]<miny) miny=fieldVerticies[1][nv];
	      nv++;
	    }
	}
    }

  /*work out the angle relative to the midpoint*/
  midx=0.5*(minx+maxx);
  midy=0.5*(miny+maxy);

  for(i=0;i<nv;i++)
    {
      fieldVerticies[0][i]-=midx;
      fieldVerticies[1][i]-=midy;
      /*args may be wrong way round*/
      fieldVerticies[2][i]=atan2(fieldVerticies[1][i],fieldVerticies[0][i]); 
    }

  /*We have the vertices of the survey polygon*/
  /*First sort them by angle*/

  do
    {
      swapped=0;
      for(i=0;i<nv-1;i++)
	{
	  if(fieldVerticies[2][i]>fieldVerticies[2][i+1])
	    {
	      t[0]=fieldVerticies[0][i]; t[1]=fieldVerticies[1][i];
	      t[2]=fieldVerticies[2][i];
	      fieldVerticies[0][i]=fieldVerticies[0][i+1];
	      fieldVerticies[1][i]=fieldVerticies[1][i+1];
	      fieldVerticies[2][i]=fieldVerticies[2][i+1];
	      fieldVerticies[0][i+1]=t[0]; fieldVerticies[1][i+1]=t[1];
	      fieldVerticies[2][i+1]=t[2];
	      swapped=1;
	    }
	}
    }
  while(swapped);

    /*Now in angle order we can work out the area*/
    fieldVerticies[0][nv]=fieldVerticies[0][0];
    fieldVerticies[1][nv]=fieldVerticies[1][0];
    fieldVerticies[2][nv]=fieldVerticies[2][0]; /*link up the ends*/

    /*calculate the area*/
    area=0.0;
    for(i=0;i<nv;i++)
      {
	eq2gal(fieldVerticies[0][i]+midx,fieldVerticies[1][i]+midy,'e',&l,&b);
	/*printf("ext %d %f %f %f %f %f\n",
	  i,fieldVerticies[0][i],fieldVerticies[1][i],
	  l*TO_DEG,b*TO_DEG,fieldVerticies[2][i]);*/
	area+=0.5*(fieldVerticies[0][i]*fieldVerticies[1][i+1]
		   -fieldVerticies[0][i+1]*fieldVerticies[1][i]);
	/*printf("%f %f\n",0.5*(fieldVerticies[0][i]*fieldVerticies[1][i+1]
	  -fieldVerticies[0][i+1]*fieldVerticies[1][i]),area);*/
      }

    area=area*TO_DEG*TO_DEG;
    /*printf("area = %f\n",area);*/

    //fprintf(logfile_ptr,"Observatory: %d Survey_Area: %f\n",obsidx,area);
  
}

double max(double x, double y)
{
  return (x>y?x:y);
};

void setupImage(struct obsfilekeywords World[], struct filekeywords *Paramfile, long* idum)
{
  int obsidx;
  //double flux_sat;

  char detfname[1000];

  
  for(obsidx = 0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      strcpy(detfname,Paramfile->obsdir);
      strcat(detfname,World[obsidx].detector);
      if(World[obsidx].im.load_detector(detfname)<0) exit(1);

      World[obsidx].im.aper.show_aperture();

      //setup the backgrounds
      World[obsidx].constbackground = pow(10,-0.4*(World[obsidx].im.background-20));
      //World[obsidx].constbackground = 0;

      World[obsidx].im.pass_seed(idum);
      if(obsidx==0)
	{   
	  if(Paramfile->outputImages && Paramfile->prettypic)
	    World[obsidx].im.set_image_properties(Paramfile->prettypicDimX, 
						  Paramfile->prettypicDimY);
	  else 
	    {
	      World[obsidx].im.minimal_image();
	    }
	}
      else
	{
	  if(abs(World[obsidx].im.psf.pixscale-World[0].im.psf.pixscale)<1e-10
	     && World[obsidx].im.aper.Naper==World[0].im.aper.Naper)
	    {
	      World[obsidx].im.set_image_properties(World[0].im.Xpix,
						    World[0].im.Ypix);
	    }
	  else
	    {
	      World[obsidx].im.set_image_properties( 
                int(ceil(World[0].im.Xpix*World[0].im.psf.pixscale / 
			 World[obsidx].im.psf.pixscale)), 
		int(ceil(World[0].im.Ypix*World[0].im.psf.pixscale /
			 World[obsidx].im.psf.pixscale))); 
	    }
	}
    }

}

void setupOrbit(struct obsfilekeywords World[], struct filekeywords *Paramfile)
{
  string ocode;
  string findstr("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");

  int elread=0; //have the elements at sim zerotime been read?

  int orbitcode;

  vector<double> data;
  vector<double> elements;
  string line;

  for(int obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
    {
      ocode = string(World[obsidx].orbitcode);

      if(ocode.find_first_of(findstr)==string::npos)
	{

	  orbitcode = atoi(World[obsidx].orbitcode);
	  //The code is a code number not a filename
	  switch(orbitcode)
	    {

	      //case0 = earth

	    case 1: //Geosynchronous
	      if(Paramfile->verbosity>0)
		cout << "Using Earth geosynch i=28 orbit for observatory " << obsidx << endl;
	      World[obsidx].orbit.resize(3);
	      World[obsidx].orbit[0].earthmoonbary();
	      World[obsidx].orbit[1].earth();
	      World[obsidx].orbit[2].geosynch(28);
	      break;

	    case 2: //L2
	      if(Paramfile->verbosity>0)
		cout << "Using L2 Lissajous orbit for observatory " << obsidx << endl;
	      World[obsidx].orbit.resize(2);
	      World[obsidx].orbit[0].earthl2();
	      //World[obsidx].orbit[1].lissajousxy(); 
	      World[obsidx].orbit[1].lissajousz(); 
	      break;

	    case 3: //Mars
	      if(Paramfile->verbosity>0)
		cout << "Using Mars orbit for observatory " << obsidx << endl;
	      World[obsidx].orbit.resize(1);
	      World[obsidx].orbit[0].mars();
	      break;

	    case 4: //Jupiter
	      if(Paramfile->verbosity>0)
		cout << "Using Jupiter orbit for observatory " << obsidx << endl;
	      World[obsidx].orbit.resize(1);
	      World[obsidx].orbit[0].jupiter();
	      break;

	    case 5: //JWST
	      if(Paramfile->verbosity>0)
		cout << "Using L2 JWST orbit for observatory " << obsidx << endl;
	      World[obsidx].orbit.resize(2);
	      World[obsidx].orbit[0].earthl2();
	      //XXX CHECK WITH MATTHEW
	      World[obsidx].orbit[1].jwst(); 
	      break;

	    case 0: 
	    default: //default to Earth's orbit
	      //check if the orbit code is less than zero
	      //if it is, the value*(-1) will be used as the phase for a jwst orbit
	      if (orbitcode<0)
		{
		  if(Paramfile->verbosity>0)
		    cout << "Using L2 JWST orbit for observatory " << obsidx << endl;
		  World[obsidx].orbit.resize(2);
		  World[obsidx].orbit[0].earthl2();
		  //XXX CHECK WIT MATTHEW
		  World[obsidx].orbit[1].jwst(-1*orbitcode/360.);

		}
	      else
		{
		if(Paramfile->verbosity>0)
		  cout << "Using Earth orbit for observatory " << obsidx << endl;
		World[obsidx].orbit.resize(2);
		World[obsidx].orbit[0].earthmoonbary();
		World[obsidx].orbit[1].earth();
		}
	    } //end switch
	}
      else
	{
	  //A filename has been passed - assume elements are in two-line form

	  string tmpf1=string(World[obsidx].orbitcode);
	  string tmpfname=string(Paramfile->obsdir) + tmpf1.substr(tmpf1.find_first_of(string("./abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")));

	  ifstream infile(tmpfname.c_str());
	  if(!infile)
	    {
	      cerr << "Error: Could not open orbit file (" << tmpfname << ") (" << World[obsidx].orbitcode << ") for observatory " << obsidx << endl;
	      exit(1);
	    }

	  while(!infile.eof())
	    {
	      getline(infile,line);
	      if(line.find_first_of(ignore)==line.npos)
		{
		  split(line,data);

		  if(data.size()==6)
		    {
		      if(elread==0)
			{
			  elements = data;
			  elread=1;
			}
		      else
			{
			  orbitalElements dummy(elements,data);
			  World[obsidx].orbit.push_back(dummy);
			  elread=0;
			} //end if elread
		    }
		  else
		    {
		      cerr << "Warning: Ignoring the following line in an orbital elements file (" << World[obsidx].orbitcode << "). Make sure there are two lines per orbit, the first with six elements, the second with their time derivatives.\nIgnoring: " << line << endl;
		    } //end if data.size=6
		} //end if ignore
	    } //end while infile

	  if(elread)
	    {
	      cerr << "Warning: There was an odd number of lines in the orbital elements file (" << World[obsidx].orbitcode << "). Make sure there are two lines per orbit, the first with six elements, the second with their time derivatives.\n" << line << endl;
	    } //end if elread
	} //end if ocode=string
    } //end for observatory
} //end function


