#ifndef FINITESRC
#define FINITESRC

#include<cmath>
#include<complex>
#include<set>
#include<stack>
#include<list>
#include<fstream>
#include<iostream>
#include<map>

#include "fsData.h"
#include "integerPowers.h"
#include "lens_base.h"
#include "src_base.h"
#include "cd.h"
#include "dcdw.h"

#define FS_DEBUG 0

///////////////////////////////////////////////////////////////////////
//
//     finiteSource
//
//        contains all the information on a single lightcurve data 
//        point
//
///////////////////////////////////////////////////////////////////////

class finiteSource
{

 private:

  map<cd,pixelContents,cdless> p; //holds rays that have been shot
  multimap<double,cd> h;   //holds values of the heuristics which decide which
                           //rays to store and discard

  int calculatefs;

 public:

  //public variables

  double t;                 //time of the data point
  cd zsc;                   //source centre
  int ptype;                //amplification used 0=single, 1=gph, 2=fs
  double allA[3];           //amplifications of each type
  lens_base* lens;            //contains information about the lens
  src_base* src;            //contains information about the source
  gphsource gphsrc;         //gph source info
  double pixSize;           //the length of a pixel side
  double pixArea;           //the area covered by a pixel

  double gphthresh;         //validity threshold for the gph approximation
  double fsthresh;          //approximate error to caluclate fs mag to

  //diagnostics

  int raysShot;
  int raysSaved;

  long pixelMemMax;         //maximum amount of memory assigned to the
                            //storage of rays - measured in units of
                            //sizeof(cd)
  long pixelMemUsed;        //current amount of memory being used to
                            //store pixels

  double hAge;              //the age of a ray

  bool drawPixelsFlag;
  bool debug;
  bool quiet;
  bool storePixels;

  int error;

  //  cd lastzsc;               //The last source position - used for making sure
                            //only short steps are taken along a track
  //  int nCalcPoints;          //Number fs points calculated so far 
                            //(used as a ray saving heuristic)

  //structures

  //functions

  finiteSource(){};

  finiteSource(src_base* src_, lens_base* lens_, double gphthresh_=1.0e-5, double fsthresh_=1.0e-3)
    {
      error = 0;
      gphthresh = gphthresh_;
      fsthresh = fsthresh_;
      src=src_;
      lens=lens_;
      gphsrc=gphsource(lens);

      debug=false;
      quiet=true;
      storePixels=false;  //do not store shot rays in memory by default

      calculatefs=2;  //default is to calculate full finite source lightcurve

      resetRayStorage();      
    };

  ~finiteSource(){};

  int calcfiniteSource(double t_);
  int calcfiniteSource(cd zsc_);
  int calcfiniteSource(int drawPix=0);
  int calcgphpoint();
  int calcgphipoint(int i);
  int magnification(cd zs, double& A, int& imgs, vector<cd>& zi, vector<double>& Ai, vector<int>& tfi);
  cd checkImage(cd image, cd zs, double& Ai, int& tfi);
  void fsFind(list<pixel>* impix);
  double fsMagnification(list<pixel>* success);
  int calcFS();
  void drawPixels(list<pixel>* pixelList, const char fname[]="pixels.txt", const char srcimgname[]="srcimg.txt");
  void drawPixels(list<pixel>* pixelList, ofstream& out, const char srcimgname[]="srcimg.txt");

  //inline functions

  inline double fsmag(cd zs, int& flag)
  {
    //returns the magnification for a given source position.
    //wrapper for more complicated functions available if the user wishes.

    flag=calcfiniteSource(zs);
    return allA[ptype];
  }

  inline void reset(double ps=5)
  {
    if(debug&&!quiet) cout << "reset()" << endl;
    pixSize=src->rs/ps;
    pixArea=sqr(pixSize);
    resetRayStorage();
    //calculateCaustic();
    error = 0;
    if(debug&&!quiet) cout << "reset() done" << endl;
  }

  inline pixel pix(cd zs)
  {
    cd zsp=zs+lens->shift;
#if FS_DEBUG
    if(debug&&!quiet) cout << "pix calculation" << endl;
#endif
    return pixel(long(floor(real(zsp)/pixSize)),long(floor(imag(zsp)/pixSize)));
  }

  inline cd xip(pixel px)
  {
    //return polar((double(floor(px/Naz))+0.5)*prad,
    //                paz*(double(px%Naz)+0.5))-lens->shift;
    return cd((double(px.x)+0.5)*pixSize,(double(px.y)+0.5)*pixSize)-lens->shift;
  }

  inline cd xip_bl(pixel px)
  {
    //return polar((double(floor(px/Naz))+0.5)*prad,
    //                paz*(double(px%Naz)+0.5))-lens->shift;
    return cd(double(px.x)*pixSize,double(px.y)*pixSize)-lens->shift;
  }

  inline bool shoot(cd zi)
  {
    double zp=abs(lens->lensEquation(zi)-zsc);
    //cout << "zp=" << zp << endl;
    if(zp<src->rs) return true;
    else return false;
  }

  inline bool shoot(cd& zi, cd& zs, double& ld)
  {
    zs = lens->lensEquation(zi);
    cd zssc=zs-zsc;
    //cout << "zp=" << zp << endl;
    if(src->inSource(zssc)) 
      {
	ld = src->limbDarkening(zssc);
	return true;
      }
    else 
      {
	ld=0.0;
	return false;
      }
  }

  inline void addNeighbours(pixel px, set<pixel>* all, list<pixel>* pending)
  {
    //add neigbouring pixels to the set of coarse pixels

#if FS_DEBUG
    if(debug&&!quiet) cout << "addNeighbours(px=" << px << "all,pending)" << endl;
#endif

    //cout << "Add neighbours: " << px << endl;
    
    pixel offsets[8]={pixel(-1,-1),pixel(0,-1),pixel(1,-1),pixel(1,0),pixel(1,1),pixel(0,1),pixel(-1,1),pixel(-1,0)};
    
    int base=0;

    pixel tpx;  
    
    for(int j=base;j<8;j++) //all other pixels
      {
	tpx=px+offsets[j];
	if(!all->count(tpx))
	  {
	    all->insert(tpx);
	    pending->push_back(tpx);
	  }
      }
  }


  //debug & miscellaneous

  void debugging() {debug=true;}
  void notDebugging() {debug=false;}
  void shh() {quiet=true;}
  void speakUp() {quiet=false;}
  void verbose() {quiet=false;}
  void saveRays() {storePixels=true;}
  void discardRays() {storePixels=false;}

  inline int nRays(int nPixels)
  {
    static const double rmsCalibration=0.3;
    //the number of rays to a side of a pixel
    //    return int(sqrt(pixArea/(src->area))/sqrt(fsthresh));
    return int(ceil( (rmsCalibration/fsthresh) * sqrt(src->area/(nPixels*pixArea)) ));
  }

  //size of memory to be used for saving pixels
  void memoryForRays(double mb) 
  {
    if(mb>0.0)
      {
	pixelMemMax=int(mb*1.0e6/sizeof(cd));
	if(debug||!quiet) cout << "Saving a maximum of " << pixelMemMax << " pixels" << endl;
      }
    else
      {
	pixelMemMax=0;
	storePixels=false;
      }
  }

  double heuristic(cd z)
  {
    return ++hAge;
  }

  void resetRayStorage()
  {
    //reset the ray storage and heuristic maps
    p.clear();
    h.clear();
    hAge=0;
    pixelMemUsed=0;
  }

  //set the type of magnification you want to calculate. Each level includes 
  //those below it
  void setfslevel(int cfs)
  {
    if(cfs<=3 && cfs>=0)
      {
	calculatefs=cfs;
      }
    else
      {
	cerr << "finiteSource: setfslevel: Error: Bad value passed to setfslevel. Choose from 0 - point source, 1 - gould approximation, 2 - full finite source with gould approximation short cut, 3 - full finite source for all points. Assuming 2." << endl;
	calculatefs=2;
      }
  }

  void pointSourceOnly() {calculatefs=0;}
  void gphOnly() {calculatefs=1;}
  void fullFiniteSource() {calculatefs=2;}

  //void calculateCaustic();

};

//////////////////////////////////////////////////////////////
//
//          structure to hold the caustics
//
//////////////////////////////////////////////////////////////

struct caustics
{
  vector<vector<cd> > caustic;
  int ncaustics;
  int topology;

  caustics(double s, double q)
  {
    caustic.clear();
    topology = causticTopology(s,q);
    ncaustics = (topology==1?3:(topology==2?1:3));
    caustic.resize(ncaustics);
  };
};

#endif /*FINITESRC*/
