#ifndef COORDS_CLASS

#include<vector>

#include "constants.h"

using namespace std;

//******************************************************************
//
//     Class for transforming between galactic and equatorial
//     coordinate systems for positions and proper motions.
//     Also converts Galactic U,V,W velocities to radial velocity
//     and proper motion based on Johnson & Soderblom 1987
//   
//     X,Y,Z  are centered on the Sun, with the GC at (+R0,0,0)
//     U +ve towards gc. V +ve towards +ve Y. W +ve towards +ve Z
//******************************************************************

class coords
{
 public:
  //transformation matrix to convert alpha,dec(J2000.0) to lII,bII
  double T[3][3];


  //transformation matrix to convert lII,bII to alpha,dec(J2000.0)
  double iT[3][3];

  void ad2lb(double a, double d, vector<double>* lb);
  void ad2lb(double a, double d, double* l, double* b);
  void lb2ad(double l, double b, vector<double>* ad);
  void lb2ad(double l, double b, double* a, double* d);

  void XYZ2Rlb(vector<double> XYZ,vector<double>* Rlb);
  void XYZ2Rlb(double X, double Y, double Z, double* R, double* l, double* b);
  void Rlb2XYZ(vector<double> Rlb, vector<double>* XYZ);
  void Rlb2XYZ(double R, double l, double b, double* X, double* Y, double* Z);

  void uvw2rvmuadlb(double l, double b, double dist, double U, double V, double W, vector<double>* rvmuadlb);
  
  double rad2hrs(double x){return x*r2h;};
  double hrs2rad(double x){return x*h2r;};
  double rad2deg(double x){return x*r2d;};
  double deg2rad(double x){return x*d2r;};

  void muad2lb(double a, double d, double mua, double mud, double* mul, double* mub);
  void mulb2ad(double l, double b, double mul, double mub, double* mua, double* mud);

  double fold(double x, double min, double max);
  
  coords()
    {

      T[0][0]=-0.0548765273295172; T[0][1]=-0.873436658110768; 
      T[0][2]=-0.483835685968435; 
      T[1][0]=0.49411066572488; T[1][1]=-0.444830355306367; 
      T[1][2]=0.746980993744103; 
      T[2][0]=-0.867665382947348; T[2][1]=-0.198076649977489; 
      T[2][2]=0.455985113757595;

      for(int i=0;i<3;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      iT[i][j]=T[j][i];
	    }
	}
    }
};

#define COORDS_CLASS
#endif
