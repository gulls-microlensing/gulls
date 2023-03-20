#ifndef FSDATA
#define FSDATA

#include "constants.h"
#include "lens_base.h"
#include "cd.h"

using namespace std;

#ifndef ULONGTYPECAST
#define ULONGTYPECAST
typedef unsigned long ulong;
#endif /*ULONGTYPECAST*/


///////////////////////////////////////////////////////////////////////
//
//     constants
//
//
///////////////////////////////////////////////////////////////////////

//directions

static const int mc=0;     //this poxel
static const int tl=1;     //middle refers to row, i.e. top, middle, bottom
static const int tc=2;     //centre refers to column, i.e. left, centre, right
static const int tr=3; 
static const int ml=4; 
static const int mr=5; 
static const int bl=6; 
static const int bc=7; 
static const int br=8; 

//gph offsets
//order centre, half plus, full plus, full cross
static const double gphOffsetX[13]={0.0, 0.5,0.0,-0.5,0.0, 1.0,0.0,-1.0,0.0, oosqrt2,-oosqrt2,-oosqrt2,oosqrt2};
static const double gphOffsetY[13]={0.0, 0.0,0.5,0.0,-0.5, 0.0,1.0,0.0,-1.0, oosqrt2,oosqrt2,-oosqrt2,-oosqrt2};

//error flags
static const int wrongNumberOfImages=1;

//thresholds
static const double magEps=1.0e-7;

//fs constants

///////////////////////////////////////////////////////////////////////
//
//     gphpoint
//
//        contains the information for a single test point of the
//        gph approximation
//
///////////////////////////////////////////////////////////////////////

struct gphpoint
{
  cd zsi;          //position of point on source plane
  int imgs;        //number of images
  double A;        //total magnification of gph point
  
  vector<cd> zi;          //positions of images on source plane
  vector<double> Ai;      //magnification of images
  vector<int> tfi;       //Are the images real? true=yes
  
  gphpoint(){};
  
  gphpoint(lens_base* lens)
  {
    //make sure there is space for each of the images
    while(int(zi.size())<lens->nImagesMax) zi.push_back(cd(0.0,0.0));
    while(int(Ai.size())<lens->nImagesMax) Ai.push_back(0.0);
    while(int(tfi.size())<lens->nImagesMax) tfi.push_back(0);
  }

};
  
///////////////////////////////////////////////////////////////////////
//
//     gphsource
//
//        contains all the information necessary for a single data
//        point of the gph amplification
//
///////////////////////////////////////////////////////////////////////
  
struct gphsource
{
  gphpoint gphpoints[13];   //data on each of the test points in the
  //gph approximation
  double Agph;              //gph magnification
  double Azeroth;           //zeroth order term in the magnification
  double Afirst;            //first    "    "   "   "        "
  double Asecond;           //second   "    "   "   "        "

  gphsource(){};

  gphsource(lens_base* lens)
  {
    //make sure the image storage is big enough
    for(int i=0;i<13;i++) gphpoints[i] = gphpoint(lens);
  }

};

///////////////////////////////////////////////////////////////////////
//
//     poxel
//
//        contains the information of a polar coordinate grid pixel,
//        or `poxel'
//
///////////////////////////////////////////////////////////////////////

struct poxel
{
  cd zi;                    //image plane ray position
  cd zs;                    //source plane landing site
  bool landed;              //whether the ray landed on the source
  bool neighbours[8];       //whether its neighbours have been tested
  int nNeighbours;          //number of neighbours checked
  int ntNeighbours;         //number of neighbours with successful rays
                            //this may not be needed
};

class pixel
{
 public:
  long x,y;

  pixel(){x=y=0;};
  pixel(long x_,long y_){x=x_; y=y_;};
  ~pixel(){};

  //inline functions

  inline bool operator<(const pixel& rhs) const
  {
    if(x<rhs.x) return true;
    else if(x>rhs.x) return false;
    else if(y<rhs.y) return true;
    else return false;
  }

  inline bool operator==(const pixel& rhs) const
  {
    return ((x==rhs.x)&&(y==rhs.y));
  }

  inline pixel& operator=(const pixel& rhs)
  {
    if(this!=&rhs)
      {
	this->x=rhs.x;
	this->y=rhs.y;
	
	return *this;
      }

    return *this;
  }

  inline pixel& operator+=(const pixel& rhs)
  {
    x+=rhs.x;
    y+=rhs.y;

    return *this;
  }

  inline pixel operator+(const pixel& other) const
  {
    pixel result=*this;
    result+=other;
    return result;
  }

  inline pixel& operator-=(const pixel& rhs)
  {
    x-=rhs.x;
    y-=rhs.y;

    return *this;
  }

  inline pixel operator-(const pixel& other) const
  {
    pixel result=*this;
    result-=other;
    return result;
  }

  inline friend ostream& operator<<(ostream& stream, pixel op)
  {
    stream << '(' << op.x << ',' << op.y << ')';
    return stream;
  }

};

//////////////////////////////////////////////////////////////
//
//
//      structure to hold the contents of a ray shot pixel
//
//
//////////////////////////////////////////////////////////////

struct pixelContents
{
  vector<cd> zs; //only need zs as z will be implied by its key in the map

  pixelContents(int nsub)
  {
    zs.resize(nsub,cd(0.0,0.0));
  }

  ~pixelContents()
  {
  }
};

#endif /*FSDATA*/
