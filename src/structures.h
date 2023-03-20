#include<vector>
#include<string>

#include "definitions.h"
#include "image.h"
#include "columnCodes.h"
#include "zodiacalLight.h"
#include "ephem.h"
#include "parallax.h"
#include "starlist.h"

#ifndef OBSSEQUENCESTRUCT
struct obssequence{
  int field;
  int nstack;
  double texp;
};
#define OBSSEQUENCESTRUCT
#endif

#ifndef OBSTRYSTRUCT
struct obsfilekeywords{
  double latitude;
  double longitude;
  double altitude;
  int finalidxnights;
  int nfields; /*NUMBER OF FIELDS FROM .CONFIG FILE */
  double readohead;
  int nepochs;
  char name[100];
  int npixx,npixy;
  double pixelsize;
  double collectingarea,primary,blockage;
  int space;
  double Aseen;   /*Area of sky over which events are seen by observatory*/
  double Aoccured; /*Area of sky over which events occured*/
  int Nseen; /*Number of events seen by observatory*/
  int Noccured; /*Number of events over entire simulated area*/
  int filter;
  int sequence_length;
  //int same_sequence;
  image im; //make sure observatory with the largest required image comes first
  image ref; //reference image for difference imaging
  double mintexp;
  int photcode; //type of photometry used

  //parameters for constructing the reference image
  double reftexp;
  double refnstack;

  //backgrounds
  zodiacalLight zodi;
  double constbackground;
  double skybackground;
  double extcoeff;

  //orbit
  vector<orbitalElements> orbit;

  struct obssequence* sequence;

  vector<double> epoch;
  //double* weatherSequence;
  vector<double> weatherSequence;
  vector<int> field;
  vector<int> weather;
  vector<double> exptime;
  vector<int> nstack;
  vector<double> jd;
  vector<double> lambda;  //ecliptic coordinate lambda of field centre
  vector<double> beta; //ecliptic coordinate beta of field centre
  vector<double> lambdasun; //ecliptic coordinate lambda the sun
  vector<double> zodiflux; //zodiacal background flux in units of mag20 per sq "

  char* weatherProfile;
  char* fieldCentreFile;
  char* observationSequence;
  char* detector;
  char* throughput;
  char* orbitcode; //codification of the orbit or filename


  //double sunriseset[NUM_SIM_DAYS+4][2];
  vector<vector<double> > sunriseset;
  /*  l,b field centres from <obsname>.config */
  double fieldCentres[2][MAX_NUM_FIELDS];
  /*  l,b field verticies computed by computeFieldVerticies */
  double fieldVerticies[2][4][MAX_NUM_FIELDS];


  obsfilekeywords()
  {
    sequence = new struct obssequence[MAX_SEQUENCE_LENGTH];
    //weatherSequence = new double[4*(NUM_SIM_DAYS+1)];

    weatherProfile = new char[1000];
    fieldCentreFile = new char[1000];
    observationSequence = new char[1000];
    detector = new char[1000];
    throughput = new char[1000];
    orbitcode = new char[1000]; 
  };

  ~obsfilekeywords()
  {
    delete[] sequence;
    //delete[] weatherSequence;

    delete[] weatherProfile;
    delete[] fieldCentreFile;
    delete[] observationSequence;
    delete[] detector;
    delete[] throughput;
    delete[] orbitcode;
  };

};
#define OBSTRYSTRUCT
#endif

#ifndef PARAMFILESTRUCT
struct filekeywords{
  
  char setseedtoclock[2];
  long Seed;
  int setseedtoclockBIT;
  double simulation_zerotime;
  int NUM_SIM_DAYS;
  int numobservatories;
  char run_name[100];
  char niterSTR[10];
  long niter;
  int principle_observatory;
  char principle_observatorySTR[10];
  char outputLightcurveSTR[10];
  double outputLightcurve; //probability of outputting a lightcurve
  int sourcecolours;
  int lenscolours;
  int Nfilters;
  double Amin;
  double large_psf_mag;
  int outputImages;
  int prettypic;
  int prettypicDimX,prettypicDimY;
  double minChiSquared; //minimum value of chi^2 for event to be interesting
  int outputOnErr;
  int outputOnDet;
  int outputOnAll;
  vector<int> validFields;
  int lenslight;
  double pllxMultiplyer; //Adjust the strength of parallax - will usually choose 0 (off) or 1 (normal)
  int verbosity;
  int choosefield;
  int identicalSequence; //reuse the lightcurve calculations with different bands
  double u0max;
  
  // A switch for appending uniform error scalings to find DeltaChi2 and n3sig cuts
  int error_scaling;

  long* seed;

  double alltime;
  double lctime;
  double phottime;
  
  char *pathdir;
  char *pathfile;
  char *basedir;
  /*char *scriptdir;
  char *paramdir;
  char *srcdir;
  char *paramfiledir;
  char *observdic;
  char *weatherdir;*/
  char *obsdir;
  char *obslist;
  char *weatherprofiledir;
  char *random_seed;
  char *simulation_zerotimeSTR;
  char *starfielddir; 
  char *starfieldlist; 
  char *sourcedir; 
  char *sourcelist; 
  char *lensdir; 
  char *lenslist; 
  char *planetdir;
  char *planetroot;
  char *outputdir;
  char *obsgroupstr;

  filekeywords()
  {
    pathdir = new char[1000];
    pathfile = new char[1000];
    basedir= new char[1000];
    obsdir = new char[1000];
    obslist = new char[1000];
    weatherprofiledir = new char[1000];
    random_seed = new char[1000];
    simulation_zerotimeSTR = new char[1000];
    starfielddir = new char[1000];
    starfieldlist = new char[1000];
    sourcedir = new char[1000];
    sourcelist = new char[1000];
    lensdir = new char[1000];
    lenslist = new char[1000];
    planetdir = new char[1000];
    planetroot = new char[100];
    outputdir = new char[1000];
    obsgroupstr = new char[1000];
  };
  
  ~filekeywords()
  {
    delete[] pathdir;
    delete[] obsdir;
    delete[] obslist;
    delete[] weatherprofiledir;
    delete[] random_seed;
    delete[] simulation_zerotimeSTR;
    delete[] starfielddir;
    delete[] starfieldlist;
    delete[] sourcedir; 
    delete[] sourcelist; 
    delete[] lensdir; 
    delete[] lenslist; 
    delete[] planetdir;
    delete[] planetroot;
    delete[] outputdir;
    delete[] obsgroupstr;
  };
};
#define PARAMFILESTRUCT
#endif

#ifndef MAGSTRUCT
class magnitudes
{
 private:

 public:
  vector<double> m;

  magnitudes(int n=0)
    {
      m.resize(n);
    }
};
#define MAGSTRUCT
#endif

#ifndef FITTED
struct fittedparams{
  double umin;
  double t0;
  double tE;
  double Fu[MAX_NUM_OBSERVATORIES];
  double Fl[MAX_NUM_OBSERVATORIES];
  double chisqvec[MAX_NUM_OBSERVATORIES];
  double chisq;
  double rS,ld1;
  double piEN, piEE;
  vector<parallax> pllx;
};
#define FITTED
#endif

#ifndef EVENTSTRUCT
struct event{

  int source, lens;
  int field;
  int id;
  //microlensing paramters
  double u0, alpha, t0, tE_h, tE_r, rE, thE, piE, piEN, piEE, rs, murel, vt, gamma; 
  //weights
  double u0max, raww, w;
  double l, b, ra, dec; //positions
  //vector<double> smag, lmag; //source and lens magnitudes
  vector<double> params; //other, non-intrinsic parameters
  // Flag array for which observatory sees the event in which field
  int isseenby[MAX_NUM_OBSERVATORIES][MAX_NUM_FIELDS]; 
  //fraction of total flux contrib by source
  double fs[MAX_NUM_OBSERVATORIES];     
  int nepochs;                            /*TOTAL NUMBER OF EPOCHS */
  int nepochsvec[MAX_NUM_OBSERVATORIES+1];    /*CUMULATIVE SUM OF EPOCHS*/
  int numobservatories;
  vector<struct fittedparams> PSPL;
  vector<struct fittedparams> FSPL;
  vector<int> flag_needFS;
  int instance;
  int lcerror;    //lightcurve generation flag
  int fisherror;  //fisher matrix calculation error
  int deterror;   //detection criteria error flag
  int detected;   //flag denoting whether event is detected

  int currentgroup;
  vector<vector<int> > obsgroups; //groups of observatories for fisher matrix calculations
  vector<string> obsgroupoutput; //Output generated for each observation group

  int allsat;
  int allsatobs[MAX_NUM_OBSERVATORIES];
  double baselineFlux[MAX_NUM_OBSERVATORIES];
  vector<int> flatlc;
  vector<double> flatchi2;
  double Amax;  //maximum measured magnification
  double umin; //minimum value of u in the lightcurve
  int peakpoint; //the epoch number of the peak point
  int outputthis;

  //parallax
  parallax pllx[MAX_NUM_OBSERVATORIES];

  //scale factors
  double scale_factor_300;
  double scale_factor_n3sig3;
  double scale_factor_n3sig6;

  //x and y position of the source in the image 
  int xsub[MAX_NUM_OBSERVATORIES]; //in sub pixel coordinates
  int ysub[MAX_NUM_OBSERVATORIES]; 
  int xpix[MAX_NUM_OBSERVATORIES]; //in pixel coordinates
  int ypix[MAX_NUM_OBSERVATORIES];   
  
  vector<double> epoch;
  vector<int> jdepoch; //points back to the point in World[obsidx].jd that corresponds to this epoch
  vector<double> alt;
  vector<int> obsidx;
  vector<int> nstack;
  vector<double> texp;
  vector<double> moonObjDist;
  vector<double> deltaVmoon;
  vector<double> Atrue;
  vector<double> Atrueerr;
  vector<double> Aobs;
  vector<double> Aerr;
  vector<double> Afit;
  vector<bool> nosat;      /*Is point unsaturated? */
  vector<double> backmag;
  vector<double> dF;
  vector<double> xs; //source position
  vector<double> ys;
  vector<double> xl1; //lens 1 position
  vector<double> yl1;
  vector<double> xl2; //lens 2 position
  vector<double> yl2;

  vector<double> data; //generic data to be output
  vector<starlist> sl;

};
#define EVENTSTRUCT
#endif

#ifndef SL_CAT
struct slcat
{
  int field;
  int start;
  int end;
  double l;
  double b;
  double dl;
  double db;
  vector<vector<double> > data;
  vector<vector<double> > mags;
  //         J     R-H     I-H     J-H           mul     mub        Vr    UU      VV      WW      Mv  CL Typ  Teff  logg Age Mass Mbol Radius [Fe/H]  l(deg)      b(deg)   RA2000.0     DEC2000.0       Dist   x(kpc)  y(kpc)  z(kpc)  Av [alpha/Fe]

  slcat()
    {
      //      start = vector<int>(MAX_STARFIELDS,-1);
      //    end = vector<int>(MAX_STARFIELDS,-1);
      //    l = vector<double>(MAX_STARFIELDS,0);
      //    b = vector<double>(MAX_STARFIELDS,0);
      //    dl = vector<double>(MAX_STARFIELDS,0);
      //    db = vector<double>(MAX_STARFIELDS,0); 

  };

  ~slcat()
  {
  
  };
};
#define SL_CAT
#endif

#ifndef P_CAT
struct pcat
{
  vector<double> data;

  pcat(vector<double>* d, int size)
  {
    data.resize(size);
    for(int i=0;i<int(d->size());i++)
      data[i] = (*d)[i];
  };
};
#define P_CAT
#endif
