#ifndef PSF_HEADER

#include<vector>
#include<string>
#include<iostream>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>

#include "integerPowers.h"

using namespace std;

/**********************************************************************

                               image + psf

                        Written by: Matthew Penny
                 Jodrell Bank Centre for Astrophysics
                        University of Manchester

                      Copyright: Matthew Penny 2011

A C++ class for simulating imaging and photometry in dense star-fields

If you wish to use this software for a publication, please contact me
before doing too much work, as I may require to be an author.
                     

**********************************************************************/

double eval_gaussian(double x, double y, void* params);

class PSF
{
 private:

  int pointsside; //the number of elements in the psf array: the number of 
                  //pixels in the kernel * the number of possible sub pixel
                  //positions
  int origin;     //the offset from index 0 of the psf center
  int psfsize;    //number of points on each side of the grid of the psf 
                  //evaluations
  int kernside;   //Side length of the kernel in pixels
  int kernsize;   //Total number of pixels in the kernel

  int spointsside; //the number of elements in the s psf array: the number of 
                   //pixels in the kernel * the number of possible sub pixel
                   //positions
  int sorigin;     //the offset from index 0 of the small psf center
  int spsfsize;    //number of points on each side of the grid of the s psf 
                   //evaluations
  int skernside;   //Side length of the small kernel in pixels
  int skernsize;   //Total number of pixels in the small kernel

  double sMissingFlux; //fraction of flux missed by using the smaller kernel
  double nMissedPix; //number of pixels missing in the smaller kernel

  double* psfxy; //array temporarily holding sampled values of psf
  long psfx; //size of the above array
  long psfy;
  double psfh; //the step size between psf samples
  double ox, oy; //the coordinate of the bottom left point of the tabulated 
                 //grid, relative to the psf peak - needs to be set by any
                 //function that loads or calculates a tabulated psf

  vector<double> psf;  //ensemble of Nsub*Nsub integrated psfs @ subpix points
  vector<double> spsf; //small psf for use on fainter stars
  vector<double> ipsf; //single psf where center can take any subpixel value

  vector<double> ctalk; //Cross-talk kernel

  vector<gsl_interp_accel*> acc; //vector of accelerators for the row splines
  vector<gsl_spline*> spline; //vector of row splines
  

  double xint(int xmin, int xmax, int y, vector<double>* p);
  double xint(int xmin, int xmax, int y, double p[]);

  int preconvolved;

 public:
  
  bool init; //flag to check if psf is initialized
  int Nkern;      //Half kernel size: there are 2Nkern+1 pixels per kernel side
  int sNkern;     //Half small kernel size
  int Nsub;       //Number of sub pixel positions
  double pixscale;     //The sidelength of each pixel

  //constructor
  PSF(int Nkern_=7, int Nsub_=8, double pixscale_=0.3)
    {
      psfxy=NULL; //don't point to anything yet

      Nkern=Nkern_;
      sNkern=int(Nkern/2);
      Nsub=Nsub_;
      pixscale=pixscale_;

      //work out some useful quantities
      kernside = 2*Nkern+1;
      skernside = 2*sNkern+1;
      kernsize = sqr(kernside);
      skernsize = sqr(skernside);
      psfsize = kernsize*sqr(Nsub);
      spsfsize = skernsize*sqr(Nsub);
      pointsside = (kernside+2)*Nsub+1;
      spointsside = (skernside+2)*Nsub+1;
      origin = (Nkern+1)*Nsub;
      sorigin = (sNkern+1)*Nsub;
      nMissedPix = kernsize - skernsize;

      preconvolved=0;

      init = false; //the psf needs recalculating
    };
  
  ~PSF()
    {
      if(psfxy!=NULL) delete[] psfxy;
      for(int i=0;i<int(spline.size());i++) gsl_spline_free(spline[i]);
      for(int i=0;i<int(acc.size());i++) gsl_interp_accel_free(acc[i]);
    };

  //the following pixval functions are dependent on the implementation of image

  //return the value of the pixel (Px,Py) for a star at subpixel coordinates (Sx,Sy)
  double pixval(int Px, int Py, int Sx, int Sy)
  {
    return psf[kernsize*(Nsub*Sy+Sx) + kernside*(Py+Nkern) + (Px+Nkern)];
  }

  //for optimization
  double pixval(int n)
  {
    return psf[n];
  }

  //for optimization
  int psfpixidx(int Px, int Py, int Sx, int Sy)
  {
    return kernsize*(Nsub*Sy+Sx) + kernside*(Py+Nkern) + (Px+Nkern);
  }

  double spixval(int Px, int Py, int Sx, int Sy)
  {
    return spsf[skernsize*(Nsub*Sy+Sx) + skernside*(Py+sNkern) + (Px+sNkern)];
  }

  //for optimization
  double spixval(int n)
  {
    return spsf[n];
  }

  //for optimization
  int spsfpixidx(int Px, int Py, int Sx, int Sy)
  {
    return skernsize*(Nsub*Sy+Sx) + skernside*(Py+sNkern) + (Px+sNkern);
  }

  double ipixval(int Px, int Py, int Sx, int Sy)
  {
    return ipsf[kernside*(Py+Nkern) + (Px+Nkern)];
  }

  //generate a psf from a function
  int generate_psf(double (*psffunc)(double, double, void*), void* params, double step);

  //load a psf from a text file
  int load_txt(string fname, double step, double centx_=1e30, double centy_=1e30, int nx=-1, int ny=-1);

  //free input psf memory
  void flush_psfxy()
  {
    delete[] psfxy;
    psfxy = NULL;
  };

  //free interpolation memory
  void flush_splines()
  {
    for(int i=spline.size()-1;i>=0;i--)
      {
	gsl_spline_free(spline[i]);
	gsl_interp_accel_free(acc[i]);
	spline.pop_back();
	acc.pop_back();
      }
    spline.clear();
    acc.clear();
  }

  //void integrate(double (*psffunc)(double, double, void*), void* params);
  //integrate the psf over a pixel
  int integrate();
  //double integrate_x(double x0, double x1, int y);
  //double integrate_pixel(double x0, double y0, double x1, double y1);
  void interpolate();
  void integrate_kernel(double x, double y, int xs=-1, int ys=-1);
  double integrate_spline(double x0, double y0, double x1, double y1);
  void extract_kernel(double x, double y, int xs, int ys);
  double extract_spline(double x0, double y0);

  //int integrate_fits(string fname, int subsamp);

  void output_psf(int si, int sj, string filename, bool overwrite=true);

  //return the amount of flux that is missed by the small aperture
  inline double missing_flux()
  {
    return sMissingFlux;
  }

  int write_psf(string filename);
  int read_psf(string filename);

  int load_crosstalk(string ifname);
  int load_crosstalk(vector<double> ct);
  void apply_crosstalk();

  void print_vars()
  {
       cout << "spsf vars: " << spointsside << " " << sorigin << " " << spsfsize << " " << skernside << " " << skernsize << " " << sMissingFlux << " " << nMissedPix << endl;
      cout << "psf vars: " << pointsside << " " << origin << " " << psfsize << " " << kernside << " " << kernsize << endl;
  }
  
};

#define PSF_HEADER
#endif
