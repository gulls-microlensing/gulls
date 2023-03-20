#include "gulls.h"
#include "structures.h"


char input_filename[100];        /* FILE HANDLING */
char output_filename[100];
char log_filename[100];
char eventprefix[100];
char instance[100];
FILE *infile_ptr;                  
FILE *outfile_ptr;
FILE *logfile_ptr;

char str[100];

struct filekeywords Paramfile;
struct obsfilekeywords World[MAX_NUM_OBSERVATORIES];
struct event Event;
struct hist BlendData[3][3]; 
struct galaxy Galaxy;

int idx, errval;
time_t st,ft,st1; 

long var = -1;
long *idum;
int ind;

double pixelDisc(double rad, double ps);
double sqr(double x){return x*x;}
double qAdd(double x, double y)
{
  if(x>y) return fabs(x)*sqrt(1+sqr(y/x));
  else return fabs(y)*sqrt(1+sqr(x/y));
}

static void usage(int status) {
  fprintf(stderr, "Usage: exigere  -i <infile> -s <instance>\n");
  exit(status);
}


int main(int argc, char *argv[]){                   /* BEGIN MAIN */

  if (argc < 3) usage(0);

  int debug=0;

  int option;

  long seed=1;
  idum = &seed;

  errval=0;

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

  if(debug) printf("readParamfile\n");
  readParamfile(input_filename, &Paramfile);
  sprintf(output_filename, "%s%s_%s.out", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation output filename */
  sprintf(log_filename, "%s%s_%s.log", Paramfile.outputdir, Paramfile.run_name, instance);   /* Create simulation log filename */
  sprintf(eventprefix, "%s%s_%s_lc", Paramfile.outputdir, Paramfile.run_name, instance);        /* Create lightcurve filename prefix */

  outfile_ptr = fopen(output_filename,"w");
  logfile_ptr = fopen(log_filename,"w");


  if (outfile_ptr == NULL)
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

  if (logfile_ptr == NULL)
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

  /* This function creates the world of telescopes */
  if(debug) printf("buildWorld\n");
  buildWorld(&Paramfile, World,idum,logfile_ptr);
  if(debug) printf("world built\n");

  /* For each observatory/filter combination work out its photometric accuracy*/
  /* as a function of magnitude */
  int obsidx;
  double mback; /* hard coded for the moment */
  double Ns, Nback, Ntot, Fs, Fback, pixel_Fback, Fbackpua;
  double sigma_phot,sigma_sys,sigma;
  double ms;


  for(obsidx=0;obsidx<Paramfile.numobservatories;obsidx++)
    {
      mback = World[obsidx].background;
      mag2flux(mback, Paramfile.primaryColour[0], 0.0, 0.0, 0.0, &Fbackpua);
      pixel_Fback = Fbackpua * pow(World[obsidx].pixelsize,2);
      Fback = pixel_Fback * pixelDisc(World[obsidx].meanseeing,World[obsidx].pixelsize);
      Nback = Fback * World[obsidx].collectingarea * World[obsidx].texp;

      printf("#Observatory %d:\n\n",obsidx);

      for(ms=15.0; ms<25.0; ms+=0.01)
	{
	  mag2flux(ms, Paramfile.primaryColour[0], 0.0, 0.0, 0.0, &Fs);

	  Ns = Fs * World[obsidx].collectingarea * World[obsidx].texp;

	  Ntot = Ns + Nback;

	  sigma_phot = sqrt(Ntot);
	  sigma_sys = World[obsidx].sigma_sys*Ntot; 
	  
	  sigma = sigma_phot * sqrt(1.0 + pow(sigma_sys/sigma_phot,2)) / Ns;

	  printf("%f %f %f\n",ms,sigma,1.0/sigma);

	}
    }
     
  return(0);
  
}

double pixelDisc(double rad, double ps)
{
  /*calculate the number of pixels in a disc of radius r*/

  /*pixScale (ps)is pixel width*/

  /*assume disc centered on central pixel*/

  int largestXY=ceil(rad/ps)+1;

  double l,r,t,b;
  double xc,yc;

  int nPix=0;
  int in;
  int i,j;

  for(i=-largestXY; i<=largestXY;i++)
    {
      xc = i*ps;
      l = xc-0.5*ps;
      r = xc+0.5*ps;
      for(j=-largestXY; j<=largestXY;j++)
	{
	  in=0;
	  yc=j*ps;
	  b = yc-0.5*ps;
	  t = yc+0.5*ps;
	  
	  in |= qAdd(b,l)<=rad;
	  in |= qAdd(b,r)<=rad;
	  in |= qAdd(t,l)<=rad;
	  in |= qAdd(t,r)<=rad;
	  if(in) nPix++;
	}
    }

  return nPix;

}
