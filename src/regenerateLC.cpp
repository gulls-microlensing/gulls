#include<vector>

#include "gulls.h"
#include<iostream>
#include<fstream>
#include<ctime>
#include<sys/timeb.h>

using namespace std;


char input_filename[100];        /* FILE HANDLING */
char output_filename[100];
char log_filename[100];
char eventprefix[100];
char instance[100];
FILE *infile_ptr;                  
ofstream outfile_ptr;
ofstream logfile_ptr;

char str[100];

struct filekeywords Paramfile;
struct obsfilekeywords* World;
struct event Event;
//struct galaxy Galaxy;
struct slcat Sources;
struct slcat Lenses;
vector<struct pcat> Planets;

//stored list of stellar magnitues
//starfield[besancon field number][level number][star number][band]
vector<vector<vector<vector<double> > > >* starfield;
vector<vector<double> >* starfieldData;

int idx;
time_t st,ft,st1; 

struct timeb mtime;

long var = -1;
long *idum;
int ind;

static void usage(int status) {
  fprintf(stderr, "Usage: exigere  -i <infile> -s <instance>\n");
  exit(status);
}


int main(int argc, char *argv[]){                   /* BEGIN MAIN */

  if (argc < 3) usage(0);

  World = new struct obsfilekeywords[MAX_NUM_OBSERVATORIES];
  starfield = new vector<vector<vector<vector<double> > > >();
  starfieldData = new vector<vector<double> >();
  /*  Event = new struct event;
  Galaxy = new struct galaxy;
  Paramfile = new struct filekeywords;*/

  int debug=0;

  int option;

  system("clear");
  printf("Exigere v0.1\n");
  printf("June 2008 \n\n");
  st=time(0); st1=st;
  printf("Execution begins: %s\n",ctime(&st));

  /* Parse command line options */

  while((option = getopt(argc, argv, "i:s:")) != EOF)
    {
      switch(option)
	{
	  
	case 'i' :
	  strcpy(input_filename,optarg);
	  break;
      
	case 's' :
	  strcpy(instance,optarg);
	  break;


	default :
	  printf("Error parsing command line arguments\n");
	  exit(1);
	}
    }


  infile_ptr = fopen(input_filename,"r");
  if (infile_ptr == NULL)
    {
      sprintf(str,"Unable to open input file: %s",input_filename);
      fmtline(str,WIDTH,"FAILED");
      exit(1);
    }
  else 
    {
      sprintf(str,"Input file: %s",input_filename);
      fmtline(str,WIDTH,"READY");
    }

  if(debug) {printf("readParamfile\n"); fflush(stdout);}
  readParamfile(input_filename, &Paramfile);
  sprintf(output_filename, "%s%s_%s.out", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation output filename */
  sprintf(log_filename, "%s%s_%s.log", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation log filename */
  sprintf(eventprefix, "%s%s_%s_lc", Paramfile.outputdir, Paramfile.run_name, instance);        /* Create lightcurve filename prefix */

  outfile_ptr.open(output_filename);
  logfile_ptr.open(log_filename);


  if (!outfile_ptr)
    {
      sprintf(str,"Unable to open output file: %s",output_filename);
      fmtline(str,2*WIDTH,"FAILED");
      exit(1);
    }
  else 
    {
      sprintf(str,"Output file: %s",output_filename);
      fmtline(str,2*WIDTH,"READY");
    }

  if (!logfile_ptr)
    {
      sprintf(str,"Unable to open output file: %s",log_filename);
      fmtline(str,2*WIDTH,"FAILED");
      exit(1);
    }
  else 
    {
      sprintf(str,"Log file: %s",log_filename);
      fmtline(str,2*WIDTH,"READY");
    }
 
 
  /* Initialise and warmup random number generator */
  idum = &var;        

  if(Paramfile.setseedtoclockBIT==1)
    {
      ftime(&mtime);
      var = -(long(mtime.time)*1000 + mtime.millitm + 14631*atoi(instance)); 
      sprintf(str,"Random seed: %ld",-var);
      fmtline(str,WIDTH,"CPU CLOCK");
    }
  else
    {
      var = Paramfile.Seed;
      sprintf(str,"Random seed: %ld",var);
      fmtline(str,WIDTH,"PARAMFILE");
    }

  ran2(idum);
  var = -var;
  gasdev(idum);

  /* Warmup random number generator */
  for (ind=0;ind<100;ind++)
    {
      ran2(idum);
      gasdev(idum);
    }

  fflush(stdout);

  /* This function creates the world of telescopes */
  if(debug) {printf("buildWorld\n"); fflush(stdout);}
  buildWorld(&Paramfile, World,idum,logfile_ptr);
  if(debug) {printf("world built\n"); fflush(stdout);}

  fflush(stdout);

  //Read in starfields
  if(debug) {printf("readStarfields\n"); fflush(stdout);}
  if(readStarfields(&Paramfile, starfield, starfieldData)==0)
    {
      sprintf(str,"Error reading starfields. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(debug) {printf("Starfields read\n"); fflush(stdout);}

  //Read in Sources
  if(debug) {printf("readSources\n"); fflush(stdout);}
  if(readSLList(true,&Paramfile, &Sources)==0)
    {
      sprintf(str,"Error reading sources. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(debug) {printf("Sources read\n"); fflush(stdout);}

  //Read in Lenses
  if(debug) {printf("readLenses\n"); fflush(stdout);}
  if(readSLList(false,&Paramfile, &Lenses)==0)
    {
      sprintf(str,"Error reading lenses. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(debug) {printf("Lenses read\n"); fflush(stdout);}

  //get a list of all available fields
  if(getValidFields(&Paramfile, &Sources, &Lenses, starfieldData)==0)
    {
      cerr << "Error: No valid fields were loaded" << endl;
      exit(1);
    }

  cout <<  Paramfile.validFields.size() << " fields read: " << Sources.end[Paramfile.validFields.back()] << " sources, " << Lenses.end[Paramfile.validFields.back()] << " lenses" << endl;

  //Read in the planets
  if(debug) {printf("readPlanets\n"); fflush(stdout);}
  if(readPlanets(&Paramfile, &Planets, instance)==0)
    {
      sprintf(str,"Error reading planets. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
      

  if(debug) {printf("about to start - for each galaxy model event\n");
    fflush(stdout);}

  /* For every Galaxy model event */
  for(idx=0; idx<int(Planets.size()); idx++)
    {
      /* Read in event parameters */
      if(debug) {printf("buildEvent\n"); fflush(stdout);}
      buildEvent(&Event, World, starfield, starfieldData,  
		 &Paramfile, &Sources, &Lenses, &Planets, idx, instance, idum); 
      if(debug) {printf("Event built\n"); fflush(stdout);}

      /* Compute the observation epochs */ 
      if(debug) {printf("timeSequencer\n"); fflush(stdout);}
      timeSequencer(World, &Event, &Paramfile, &Sources, &Lenses);  
      if(debug) {printf("time sequenced\n"); fflush(stdout);}
 
      /* Generate lightcurve */
      if(debug) {printf("lightcurveGenerator\n"); fflush(stdout);}
      lightcurveGenerator(&Paramfile, &Event, World, &Sources, &Lenses,
			  logfile_ptr);
      if(debug) {printf("lightcurve generated\n"); fflush(stdout);}

      //Did we detect what we are interested in?
      if(debug) {printf("detectionCriteria\n"); fflush(stdout);}
      detectionCuts(&Paramfile, &Event, World, &Sources, &Lenses);
      if(debug) {printf("detection criteria applied\n"); fflush(stdout);}

      //Output the lightcurve if desired
      if(debug){printf("outputLightcurve\n"); fflush(stdout);}
      outputLightcurve(&Event,&Paramfile,&Sources,&Lenses);
      if(debug){printf("lightcurve ouput\n"); fflush(stdout);}

      //Output images if desired
      if(debug){printf("output images\n"); fflush(stdout);}
      outputImages(&Event, World, &Sources, &Paramfile);
      if(debug){printf("images outputted\n"); fflush(stdout);}

      //Write out the events parameters and data to the appropriate file
      if(Event.lcerror || Event.deterror)
	{
	  if(Event.lcerror)
            sprintf(str,"\nDiscarding event %d (Failed lightcurve generation)",
		    idx);
	  if(Event.deterror)
            sprintf(str,"\nDiscarding event %d (Failed detection criteria)",
		    idx);
	  fmtline(str,WIDTH,"OKAY"); 
	  writeEventParams(&Paramfile, &Event, &Sources, &Lenses, logfile_ptr);
	}
      else //otherwise
	{
	  writeEventParams(&Paramfile, &Event, &Sources, &Lenses, outfile_ptr);
	}

    } //end for each event

  if (fclose(infile_ptr) == -1)
    {
      sprintf(str,"\nClosing files");
      fmtline(str,WIDTH,"FAILED");
      exit(1);
    }
  else 
    {
      sprintf(str,"\nClosing files");
      fmtline(str,WIDTH,"DONE");
    }

  delete[] World;
  delete starfield;
  delete starfieldData;

  st = time(0);
  printf("Execution ends: %s",ctime(&st));
  clock2str(st,st1);
  
  return(0);
  
}
