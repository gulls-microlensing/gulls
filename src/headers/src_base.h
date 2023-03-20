#ifndef SOURCE_BASE
#define SOURCE_BASE

#include "cd.h"

///////////////////////////////////////////////////////////////////////
//
//     src_base
//
//        Abstract base class that contains all the information 
//           on the source, and a function to calculate the limb 
//           darkened flux at a given source point. The constructor 
//           should calculate the gph limb darkening coefficients
//           and a source area. The limb darkening should be defined
//           such that the integral over the source face of the ld
//           profile is 1.
//        
//
///////////////////////////////////////////////////////////////////////

class src_base
{

 private:

 public:

  double rs;                //source radius
  double ldGamma;           //limb darkening parameter

  double area;              //area of the source

  double twoCoeff;          //second order coefficient for gph calculation
  double fourCoeff;         //fourth order coefficient for gph calculation

  //src_base(){}; 
  //initializes all the variables of the sourceInfo class to defaults
  //src_base(double rs_, double ldGamma_, double fb_){};
  //initializes all the variables of the sourceInfo class

  virtual double limbDarkening(cd& zssc) = 0;
  //returns the value of the limb darkening at a point zssc - measured from
  //the source centre

  virtual bool inSource(cd& zsc,cd& zs) = 0;
  virtual bool inSource(cd& zssc) = 0;
  //is the point zs inside the source positioned at zsc

  virtual ~src_base(){}; //destructor
};


#endif /*SOURCE_BASE*/
