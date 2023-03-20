#include<vector>

#include "gulls.h"
#include<iostream>
#include<fstream>
#include<ctime>
#include<sys/timeb.h>

using namespace std;


char input_filename[1000];        /* FILE HANDLING */
char output_filename[1000];
char log_filename[1000];
char eventprefix[1000];
char instance[100];
FILE *infile_ptr;                  
ofstream outfile_ptr;
ofstream logfile_ptr;

char str[1100];

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
  fprintf(stderr, "Usage: mabuls  -i <infile> -s <instance> {-f <field>} {-d}\n");
  exit(status);
}


int main(int argc, char *argv[]){                   /* BEGIN MAIN */

  //set up program timers
  Paramfile.alltime=0;
  Paramfile.lctime=0;
  Paramfile.phottime=0;

  double tbuildevent=0, ttimesequencer=0, tdetcuts=0, tio=0, tgeneration=0, phottime=0;
  timespec allstart, allend;
  timespec tstart, tend;
  long nsec;
  clock_gettime(CLOCK_REALTIME,&allstart);

  if (argc < 3) usage(0);

  World = new struct obsfilekeywords[MAX_NUM_OBSERVATORIES];
  starfield = new vector<vector<vector<vector<double> > > >();
  starfieldData = new vector<vector<double> >();
  /*  Event = new struct event;
  Galaxy = new struct galaxy;
  Paramfile = new struct filekeywords;*/

  int debug=0;
  Paramfile.verbosity=0;

  int option;

  int field=-1;

  system("clear");
  printf("MaBuLS v0.1\n");
  printf("November 2012 \n\n");
  st=time(0); st1=st;
  printf("Execution begins: %s\n",ctime(&st));

  /* Parse command line options */

  while((option = getopt(argc, argv, "i:s:f:d")) != EOF)
    {
      switch(option)
	{
	  
	case 'i' :
	  strcpy(input_filename,optarg);
	  break;
      
	case 's' :
	  strcpy(instance,optarg);
	  break;

	case 'd' :
	  debug=1;
	  Paramfile.verbosity++;
	  break;

	case 'f' :
	  field = atoi(optarg);
	  break;


	default :
	  printf("Error parsing command line arguments\n");
	  exit(1);
	}
    }

  Paramfile.choosefield=field;

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

  if(Paramfile.verbosity) {printf("readParamfile\n"); fflush(stdout);}
  readParamfile(input_filename, &Paramfile);
  if(field<0)
    {
      sprintf(output_filename, "%s%s_%s.out", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation output filename */
      sprintf(log_filename, "%s%s_%s.log", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation log filename */
    }
  else
    {
      sprintf(output_filename, "%s%s_%s_%d.out", Paramfile.outputdir, Paramfile.run_name, instance,field);   /* Create simulation output filename */
      sprintf(log_filename, "%s%s_%s_%d.log", Paramfile.outputdir, Paramfile.run_name, instance,field);   /* Create simulation log filename */
    }

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
  Paramfile.seed = &var;

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
  if(Paramfile.verbosity) {printf("buildWorld\n"); fflush(stdout);}
  clock_gettime(CLOCK_REALTIME,&tstart);
  buildWorld(&Paramfile, World,idum,logfile_ptr);
  clock_gettime(CLOCK_REALTIME,&tend);
  nsec = tend.tv_nsec - tstart.tv_nsec;
  cout << "Buildworld: " <<  double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0)) + double(nsec<0?nsec+1000000000:nsec)*1.0e-9 << " sec" << endl;
  
  if(Paramfile.verbosity) {printf("world built\n"); fflush(stdout);}

  fflush(stdout);

  clock_gettime(CLOCK_REALTIME,&tstart);

  //Read in starfields
  if(Paramfile.verbosity) {printf("readStarfields\n"); fflush(stdout);}
  if(readStarfields(field, &Paramfile, starfield, starfieldData)==0)
    {
      sprintf(str,"Error reading starfields. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(Paramfile.verbosity) {printf("Starfields read\n"); fflush(stdout);}

  //Read in Sources
  if(Paramfile.verbosity) {printf("readSources\n"); fflush(stdout);}
  if(readSLList(field,true,&Paramfile, &Sources)==0)
    {
      sprintf(str,"Error reading sources. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(Paramfile.verbosity) {printf("Sources read\n"); fflush(stdout);}

  //Read in Lenses
  if(Paramfile.verbosity) {printf("readLenses\n"); fflush(stdout);}
  if(readSLList(field,false,&Paramfile, &Lenses)==0)
    {
      sprintf(str,"Error reading lenses. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }
  if(Paramfile.verbosity) {printf("Lenses read\n"); fflush(stdout);}

  //get a list of all available fields
  if(getValidFields(&Paramfile, &Sources, &Lenses, starfieldData)==0)
    {
      cerr << "Error: No valid fields were loaded" << endl;
      exit(1);
    }

  cout <<  Paramfile.validFields.size() << " fields read: " << Sources.end[Paramfile.validFields.back()] << " sources, " << Lenses.end[Paramfile.validFields.back()] << " lenses" << endl;

  //Read in the planets
  if(Paramfile.verbosity) {printf("readPlanets\n"); fflush(stdout);}
  if(readPlanets(&Paramfile, &Planets, instance, field)==0)
    {
      sprintf(str,"Error reading planets. Exiting.");
      fmtline(str,WIDTH,"NON-FATAL");
      exit(1);
    }

  clock_gettime(CLOCK_REALTIME,&tend);
  nsec = tend.tv_nsec - tstart.tv_nsec;
  tio += double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0))
    + double(nsec<0?nsec+1000000000:nsec)*1.0e-9; 
      

  if(Paramfile.verbosity) {printf("about to start - for each galaxy model event\n");
    fflush(stdout);}

  /* For every Galaxy model event */
  for(idx=0; idx<int(Planets.size()); idx++)
    {
      /* Read in event parameters */
      if(Paramfile.verbosity) {printf("buildEvent\n"); fflush(stdout);}
      clock_gettime(CLOCK_REALTIME,&tstart);
      buildEvent(&Event, World, starfield, starfieldData,  
		 &Paramfile, &Sources, &Lenses, idx, instance, idum);
      getPlanetvals(&Event, World, &Paramfile, &Sources, &Lenses, &Planets);
      clock_gettime(CLOCK_REALTIME,&tend);
      nsec = tend.tv_nsec - tstart.tv_nsec;
      tbuildevent += double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0))
    + double(nsec<0?nsec+1000000000:nsec)*1.0e-9; 
      if(Paramfile.verbosity) {printf("Event built\n"); fflush(stdout);}

      /* Compute the observation epochs */ 
      if(Paramfile.verbosity) {printf("timeSequencer\n"); fflush(stdout);}
      clock_gettime(CLOCK_REALTIME,&tstart);
      timeSequencer(World, &Event, &Paramfile, &Sources, &Lenses); 
      clock_gettime(CLOCK_REALTIME,&tend);
      nsec = tend.tv_nsec - tstart.tv_nsec;
      ttimesequencer += double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0))
    + double(nsec<0?nsec+1000000000:nsec)*1.0e-9;  
      if(Paramfile.verbosity) {printf("time sequenced\n"); fflush(stdout);}
 
      //Did we detect what we are interested in?
      if(Paramfile.verbosity) {printf("detectionCriteria\n"); fflush(stdout);}
      clock_gettime(CLOCK_REALTIME,&tstart);
      detectionCuts(&Paramfile, &Event, World, &Sources, &Lenses);
      clock_gettime(CLOCK_REALTIME,&tend);
      nsec = tend.tv_nsec - tstart.tv_nsec;
      tdetcuts += double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0))
    + double(nsec<0?nsec+1000000000:nsec)*1.0e-9; 
      if(Paramfile.verbosity) {printf("detection criteria applied\n"); fflush(stdout);}

      //Output images if desired
      if(Paramfile.verbosity){printf("output images\n"); fflush(stdout);}
      outputImages(&Event, World, &Sources, &Paramfile);
      if(Paramfile.verbosity){printf("images outputted\n"); fflush(stdout);}

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
      clock_gettime(CLOCK_REALTIME,&tend);
      nsec = tend.tv_nsec - tstart.tv_nsec;
      tio += double((tend.tv_sec - tstart.tv_sec) - (nsec<0?1:0))
	+ double(nsec<0?nsec+1000000000:nsec)*1.0e-9; 

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

  clock_gettime(CLOCK_REALTIME,&allend);
  nsec = allend.tv_nsec - allstart.tv_nsec;
  Paramfile.alltime = double((allend.tv_sec - allstart.tv_sec) - (nsec<0?1:0))
    + double(nsec<0?nsec+1000000000:nsec)*1.0e-9;

  cout << "Program timings:\n";
  cout << "Total:          " << Paramfile.alltime << " sec\n";
  cout << "LightcurveGen:  " << tgeneration << " sec\n";
  cout << "Photometry:     " << phottime << " sec\n";
  cout << "Build event:    " << tbuildevent << " sec\n";
  cout << "Time sequencer: " << ttimesequencer << " sec\n";
  cout << "Detection cuts: " << tdetcuts << " sec\n";
  cout << "File io:        " << tio << " sec\n";

  st = time(0);
  printf("Execution ends: %s",ctime(&st));
  clock2str(st,st1);
  
  return(0);
  
}
