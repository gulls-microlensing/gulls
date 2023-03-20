#ifndef IMAGE_HEADER

#include<vector>
#include<cmath>
#include<iostream>
#include<ctime>
#include<algorithm>
#include<limits.h>

#include "constants.h"
#include "integerPowers.h"
#include "psf.h"
#include "aperture.h"
#include "starlist.h"
#include "random.h"

using namespace std;

/**********************************************************************

                               image + psf

                        Written by: Matthew Penny
                 Jodrell Bank Centre for Astrophysics
                        University of Manchester

                      Copyright: Matthew Penny 2011-2013

A C++ class for simulating imaging and photometry in dense star-fields

If you wish to use this software for a publication, please contact me
before doing too much work, as I may require to be an author.
                     

**********************************************************************/

//Image position codes

const int NUM_IMAGE_CODES = 2;
const int IMAGE_CENTER = -1;
const int IMAGE_CENTRE = -1; //duplicate code for spelling - don't count twice
const int SUBIMAGE_CENTER = -2;
const int SUBIMAGE_CENTRE = -2;

class image
{
 private:

  long* seed;  //random seed

 public:

  //image propeties
  int Xpix, Ypix;        //Number of pixels per side
  int Npix;              //Total number of pixels - N^2
 
  //detector parameters
  double bias_level;   //Number of counts per pixel in bias image
  double read;         //readnoise in counts per read per pixel
  double dark;         //dark current in counts per sec per pix
  double therm;        //thermal noise in counts per sec per pix
  double zero;         //no of counts from a reference magnitude
  double zeromag;      //magnitude of the zero point
  long depth;          //the dynamic range of the detector from 0-depth
  double diameter;     //telescope diameter
  double blockage;     //tube blockage diameter
  double systematic;   //systematic error applied to total photometry
  double pixsize;      //physical size of the pixel (in microns)
  double crflux;       //cosmic ray flux (CR m-1 s-1)
  double pixcrflux;    //cosmic ray flux per pixel
  double fullwell;     //Maximum number of electrons that can be held
  double gain;         //Inverse gain - e-/ADU
  int bleeding;        //Does charge bleed from the pixels?
  double bleedacross;  //What fraction of charge bleeds across rows

  double largepsfmag;  //magnitude brighter than which a larger psf is used

  //sky properties
  double background;   //background in mags per sq arcsec
  double psfback;      //background due to flux missing from the small psf
  double backflux;     //flux in a pixel due to background

  //data properties
  double texp;         //total exposure time
  int nstack;          //number of stacked images

  //psf properties
  double fwhm;          //PSF full width half max
  PSF psf;              //The integrated psf (see psf.cpp)
  string psffile;       //filename of file containing numerical psf
  double psfscale;      //only used temporarily

  //photometry aperture properties
  aperture aper;  //aperture mask

  //images
  vector<double> timage; //the underlying true image in floating point counts
  vector<int> counts;    //the realization of the counts

  vector<double> wosrc; //true image without the source of interest

  //fast microlensing photometry storage
  double fast_src;   //the photon rate provided by the unmagnified source
  double fast_blend; //photons from blend stars
  double fast_total; //the photon rate from everything else
  double fast_satlimit; //source magnification that will saturate a pixel
  double fast_aperpix; //number of pixels in the aperture that are read

  double refflux, reffluxerror; //reference flux used in dia photometry (must be initialized)
  double irefflux, ireffluxerror; //reference flux used in dia photometry (must be initialized)

  //WCS parameters
  double CDELT1, CDELT2, CRPIX1, CRPIX2;
  double CD1_1, CD1_2, CD2_1, CD2_2;
  double CRVAL1, CRVAL2;
  string CTYPE1, CTYPE2, CUNIT1, CUNIT2;
  double MJDOBS; 

  //constructor
  image()
    {
      seed = NULL;
      //set important, easily forgotten values to their defaults
      set_bitdepth();
      //set_pixscale();
      set_background();
      set_bias_level();
      set_zeropoint();
      set_largepsfmag();
      set_gain();
      set_bleeding(0);
      
      psfback = 0;
      refflux=irefflux=reffluxerror=ireffluxerror=0;
    };
  ~image(){};

  //utility functions

  //calculate the pixel index
  inline int idx(int i, int j)
  {
    return j*Xpix + i;
  }

  //processes a code to a pixel coordinate
  inline void xycode(int code, int* x, int* y)
  {
    switch(code)
      {
      case IMAGE_CENTER:
	*x = Xpix/2; *y = Ypix/2;
	break;

      case SUBIMAGE_CENTER:
	*x = Xpix*psf.Nsub/2; *y = Ypix*psf.Nsub/2;
	break;
	
      default:
	cerr << __FUNCTION__ 
	     << ": Warning: The image code used has not been defined" << endl;
	*x = -9999; *y = -9999;
      }    
  };

  //processes a code to a pixel index
  inline int xycode(int code)
  {
    switch(code)
      {
      case IMAGE_CENTER:
	return idx(Xpix/2,Ypix/2);
	break;
	
	//subimage has no meaning in this context
      default:
	return code;
      }    
  };

  //calculate the number of counts in a real image pixel
  inline double detector_counts(double texp_, double nstack_=1)
  {
    //counts in a pixel due to the detector + telescope
    double nc=0;
    
    for(int i=0; i<nstack_; i++)
      {
	double fread = read*gasdev(seed);
	nc += (fread<0?long(fread-0.5):long(fread+0.5)) 
	  + poisson(texp_*(dark+therm),seed);
      }
  
    return nc;
  }
  
  //calculate the number of counts per pixel in the absence of noise
  inline double ideal_detector_counts(double texp_, double nstack_=1)
  {
    //read noise has zero mean so doesn't count
    return texp_*nstack_*(dark+therm);
  }

  //calculate the uncertainty on the detector counts in a pixel
  inline double detector_error2(double texp_, double nstack_=1)
  {
    return nstack_*((read*read+0.5*gain*gain) + texp_*(dark + therm));
  }

  //calculate the number of counts in an image due to sky background
  inline double bg_counts(double texp_, double nstack_=1)
  {
    //randomization left to user
    return texp_*nstack_*backflux;
  }

  //convert magnitudes to a count flux
  inline double mag2flux(double mag)
  {
    return zero * pow(10.0,-0.4*(mag-zeromag));
  };


  //initialization functions

  //Pass a pointer to the RNG seed
  void pass_seed(long* idum)
  {
    seed = idum;
  };

  //reseed the random number generator
  void reseed(long x=0);

  //Set the WCS parameters
  void set_wcs_params(double mjd, double l=0.0, double b=0.0, double rot=0.0);

  //set the number of counts in a bias image pixel
  void set_bias_level(double bias_level_=200)
  {
    bias_level=bias_level_;
  };

  //set the absolute zero point counts
  void set_zeropoint(double flux=100.0, double mag=20.0, double diameter_=2.0/sqrt(pi), double blockage_=0.0)
  {
    //assumes all throughput is accounted for
    //zero = flux;
    //zeromag = mag;
    zero = 1.0;
    zeromag = mag + 2.5*log10(flux); //more forgiving of forgetful programmers
  };

  //set the detector noise properties
  void set_detector_noise(double read_=0, double dark_=0, double therm_=0, double systematic_=0)
  {
    read = read_;
    dark = dark_;
    therm = therm_;
    systematic = systematic_;
  }

  //set the sky background in mags per sq arcsec
  void set_background(double mag_=22.0)
  {
    background = mag_;
    backflux = mag2flux(background)*sqr(psf.pixscale);
  }

  //set the image dimesnsions
  void set_image_properties(int Xpix_, int Ypix_)
  {
    Xpix = Xpix_;
    Ypix = Ypix_;
    Npix = Xpix*Ypix;

    reset_image();
    reset_detector();
  }

  void set_bitdepth(long depth_=65536)
  {
    if(depth_>=0 && depth_<1024) //assume the input is in n_bits
      depth = long(pow(2.0,double(min(depth_,long(log2(LONG_MAX))))));
    else if(depth_<0) //assume essentially unlimited
      depth = LONG_MAX;
    else
      depth = depth_;
  }

  void set_fullwell(long fullwell_=2147483647)
  {
    if(fullwell_<0) fullwell = depth;
    else fullwell = fullwell_;
  }

  void set_gain(double gain_=1)
  {
    gain = gain_;
  }

  void set_bleeding(int bleeding_=1, double bleedacross_=0.025)
  {
    bleeding = bleeding_;
    bleedacross = bleedacross_;
  }

  //set up and empty psf - even Nsub is more accurate
  void empty_psf(double fwhm_, double pixscale_=0.3, int Nkern_=-1, int Nsub_=8)
  {
    fwhm = fwhm_;
    if(Nkern_==-1)
      {
	Nkern_ = int(ceil(4*fwhm/double(pixscale_)));
      }
    psf = PSF(Nkern_,Nsub_,pixscale_);
  }

  //set up a gaussian psf
  void gaussian_psf(double fwhm_, double pixscale_=0.3, int Nkern_=-1, int Nsub_=8)
  {
    //Nkern_ is the half kernel size

    fwhm = fwhm_;
    void* param = &fwhm;

    if(Nkern_==-1)
      {
	Nkern_ = int(ceil(4.0*fwhm/double(pixscale_)));
      }
    
    psf = PSF(Nkern_,Nsub_,pixscale_);
    psf.generate_psf(&eval_gaussian,param,max(pixscale_/Nsub_,pixscale_/5));
    psf.integrate();

  }

  //generate a gaussian psf from parameters already loaded
  void gaussian_psf()
  {
    //Nkern_ is the half kernel size

    void* param = &fwhm;
    
    psf = PSF(psf.Nkern,psf.Nsub,psf.pixscale);
    psf.generate_psf(&eval_gaussian, param,
		     max(psf.pixscale/psf.Nsub,psf.pixscale/5));
    psf.integrate();

  }

  //set up a user supplied psf - even Nsub is more accurate
  void custom_psf(double fwhm_, double (*psf_function)(double, double, void*), void* params, double samp=-1, double pixscale_=0.3, int Nkern_=-1, int Nsub_=8)
  {
    //Nkern_ is the half kernel size
    //fwhm_ is needed for characterizing when Nkern left to default

    fwhm = fwhm_;
    if(Nkern_==-1)
      {
	Nkern_ = int(ceil(4.0*fwhm/double(pixscale_)));
      }
    
    psf = PSF(Nkern_,Nsub_,pixscale_);

    if(samp<0) samp = pixscale_/max(Nsub_,5);
    psf.generate_psf(psf_function,params,samp);
    psf.integrate();

  }

  //Load a tabulated psf from a file
  void load_psf(string psffile_, double psfscale_)
  {
    psf.load_txt(psffile_, psfscale_);
    psf.integrate();
  }

  //set the photometry aperture to a standard 5x5 square with corners missing
  void set_standard_aperture()
  {
    aper.five_nocorners();
    aper.show_aperture();
  }

  void set_circular_aperture(double r)
  {
    aper.generate_aperture(r, psf.pixscale);
    //aper.show_aperture();
  }

  //set the photometry aperture to a standard 5x5 square
  void set_square_aperture(int sides=5)
  {
    aper.square(sides);
  }

  //set cosmic ray parameters
  void set_cr_parameters(double pixsize_, double crflux_)
    {
      pixsize = pixsize_;
      crflux = crflux_;

      pixcrflux = crflux*sqr(pixsize*1e-6); //convert to pixel flux
    }

  //load the detector parameters from a file
  int load_detector(string filename);

  //set up a minimally sized image for photometry
  void minimal_image(); 

  //set the magnitude below which images are added with the full psf
  void set_largepsfmag(double mag_=25.0)
  {
    largepsfmag = mag_;
  }


  //image manipulation functions

  //add a star at a random position
  inline void addstar(double mag)
  {    
    addstar(-psf.Nkern*psf.Nsub, (Xpix+psf.Nkern)*psf.Nsub, 
	    -psf.Nkern*psf.Nsub, (Ypix+psf.Nkern)*psf.Nsub, mag);
  };
  int addstar(double mag, starlist* sl);


  //add a star at a random position in the specified range
  inline void addstar(int xmin, int xmax, int ymin, int ymax, double mag)
  {
    addstar(randint(xmin,xmax,seed),randint(ymin,ymax,seed),mag);
  }

  inline int addstar(int xmin, int xmax, int ymin, int ymax, double mag, starlist* sl)
  {
    int x = randint(xmin,xmax,seed);
    int y = randint(ymin,ymax,seed);
    return addstar(x,y,mag,sl);
  }



  //add a star at a special position
  inline void addstar(int code, double mag, bool sub=false) 
  {
    int x,y;
    xycode(code,&x,&y);
    addstar(x,y,mag,sub);
  }
  inline int addstar(int code, double mag, starlist* sl)
  {
    int x,y;
    xycode(code,&x,&y);
    return addstar(x,y,mag,sl);
  }

  //add a star at a given position
  bool addstar(int x, int y, double mag, bool sub=false);
  inline int addstar(int x, int y, double mag, starlist* sl)
  {
    int added = addstar(x,y,mag);
    if(added)
      {
	sl->x.push_back(x);
	sl->y.push_back(y);
	sl->mag.push_back(mag);
	sl->nstars++;
      }
    return added;
  }

  //add star at general position
  void freeaddstar(double x, double y, double mag, bool sub=false);
  inline int freeaddstar(double x, double y, double mag, freestarlist* sl)
  {
    sl->x.push_back(x);
    sl->y.push_back(y);
    sl->mag.push_back(mag);
    sl->nstars++;
    freeaddstar(x,y,mag);
    return sl->x.size()-1;
  }

  inline void freeaddstar(double mag)
  {    
    freeaddstar(-psf.Nkern*psf.Nsub, (Xpix+psf.Nkern)*psf.Nsub, 
	    -psf.Nkern*psf.Nsub, (Ypix+psf.Nkern)*psf.Nsub, mag);
  };
  int freeaddstar(double mag, starlist* sl);


  //add a star at a random position in the specified range
  inline void freeaddstar(double xmin, double xmax, double ymin, double ymax, double mag)
  {
    freeaddstar(xmin+(xmax-xmin)*ran2(seed),ymin+(ymax-ymin)*ran2(seed),mag);
  }

  inline void freeaddstar(double xmin, double xmax, double ymin, double ymax, double mag, freestarlist* sl)
  {
    double x = xmin+(xmax-xmin)*ran2(seed);
    double y = ymin+(ymax-ymin)*ran2(seed);
    freeaddstar(x,y,mag,sl);
  }


  //subtract a star at a given position
  inline void substar(int x, int y, double mag)
  {
    addstar(x,y,mag,true);
  };

  //subtract a star at a given position
  inline void freesubstar(double x, double y, double mag)
  {
    freeaddstar(x,y,mag,true);
  };

  //subtract a star at a special position
  inline void substar(int code, double mag) 
  {
    int x,y;
    xycode(code,&x,&y);
    substar(x,y,mag);
  }


  //subtract a star from a star list, adding it to a deleted list
  int substar(int index, starlist* from, starlist* to);

  //subtract a star from a star list, losing all knowledge of it
  void substar(int index, starlist* from);

  //pop a star from the back of a list and remove from image
  inline void popstar(starlist* from)
  {
    substar(from->x.back(),from->y.back(),from->mag.back());
    from->x.pop_back();
    from->y.pop_back();
    from->mag.pop_back();
  }

  //keep hold of the star in a different list
  int popstar(starlist* from, starlist* to); 

  //add the background due to truncated psf tails
  double addpsfbg(int xmin, int xmax, int ymin, int ymax);

  //remove the background due to missed psf tails
  void subpsfbg(double Ncounts);

  //add a field of stars from a file
  int addfield(double solid_angle, string filename, int column, starlist* sl);
  inline int addfield(double solid_angle, string filename, int column)
  {
    return addfield(solid_angle, filename, column, NULL);
  }

  //add a field of stars from a list of stars
  int addfield(double solid_angle, vector<double>* mags, starlist* sl);
  inline int addfield(double solid_angle, vector<double>* mags)
  {
    return addfield(solid_angle, mags, NULL);
  }

  //add a field of stars from a luminosity function
  int addfield(vector<double>* mag, vector<double>* density, starlist* sl);
  inline int addfield(vector<double>* mag, vector<double>* density)
  {
    return addfield(mag, density, NULL);
  }

  //add a list of stars to the image
  inline void addstarlist(starlist* sl)
  {
    for(int i=0;i<int(sl->x.size());i++)
      {
	addstar(sl->x[i],sl->y[i],sl->mag[i]);
      }
  }

  //calculate the sub-pixel dimensions over which to add a star field
  //void field_dimensions(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield, int& nrepeats);
  void field_dimensions(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield);

  //Same as field dimensions, but incorporates the possiblilty of roll
  void field_dimensions_roll(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield);


  //WCS functions
  inline void set_time(double mjdobs)
  {
    MJDOBS = mjdobs;
  }

  inline void set_crpix(double x, double y)
  {
    CRPIX1=x;
    CRPIX2=y;
  }

  inline void set_crval(double l, double b)
  {
    CRVAL1=l;
    CRVAL2=b;
  }

  inline void set_wcs_rot(double rot)
  {
    CD1_1 = CDELT1*cos(rot);
    CD1_2 = -CDELT1*sin(rot);
    CD2_1 = CDELT1*sin(rot);
    CD2_2 = CDELT1*cos(rot);
  }
  
  //imaging functions

  //clear all counts from the detector
  void reset_detector(); 

  //clear the true image
  void reset_image();

  //expose the detector for texp_ seconds and stack nstack_ times
  void expose(double texp_, int nstack_=1);
  void ideal_expose(double texp_, int nstack_=1);
  void bleedcharge(vector<long> &charge);

  //write the measured image to a fits file
  int write_fits(string filename, bool overwrite=false, bool subbias=true);

  //write the true image to a fits file
  int write_truefits(string filename, bool overwrite=false);

  //create a test image with stars placed at various sub-pixel positions
  double sub_pixel_test(string filename, bool writefits=true);
  inline double sub_pixel_test()
  {
    return sub_pixel_test(string("/dev/null"),false);
  };


  //photometry functions

  //perform photometry at image position (x,y), on the currently stored image 
  void photometry(int x, int y, double* ncounts, double* error, int* satflag, aperture* ap=NULL);
  inline void photometry(int code, double* ncounts, double* error, int* satflag, aperture* ap=NULL)
  {
    int x,y;
    xycode(code,&x,&y);
    photometry(x, y, ncounts, error, satflag, ap);
  }

  //perform photometry at image position (x,y), on a stack of nstack_ images 
  //of exposure time texp_
  //optomized for photometry of a single star, but not suitable for photometry
  //of multiple stars in a single image
  void quick_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, aperture* ap=NULL);
  inline void quick_photometry(int code, double texp_, int nstack_, double* ncounts, double* error, int* satflag, aperture* ap=NULL)
  {
    int x,y;
    xycode(code,&x,&y);
    quick_photometry(x, y, texp_, nstack_, ncounts, error, satflag, ap);
  }

  //perform photometry with no scatter at image position (x,y), on a stack of
  //nstack_ images of exposure time texp_ - use backg=false if the image has no
  //added background (i.e. if addbg() has not been called)
  void ideal_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg=true, aperture* ap=NULL);
  inline void ideal_photometry(int code, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg=true, aperture* ap=NULL)
  {
    int x,y;
    xycode(code,&x,&y);
    ideal_photometry(x, y, texp_, nstack_, ncounts, error, 
		     satflag, backg, ap);
  }

  /*void weighted_ideal_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg, aperture* ap);
   inline void weighted_ideal_photometry(int code, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg=true, aperture* ap=NULL)
  {
    int x,y;
    xycode(code,&x,&y);
    ideal_photometry(x, y, texp_, nstack_, ncounts, error, 
		     satflag, backg, ap);
		     }*/

  //perform photometry at image position (x,y), on a stack of nstack_ images 
  //of exposure time texp_
  //optomized for photometry of a single star, but not suitable for photometry
  //of multiple stars in a single image
  //Calculates both ideal and scattered photometry
  void ideal_and_scattered_photometry(int x, int y, double texp_, int nstack_, double* ncounts_i, double* error_i, double* ncounts_s, double* error_s, int* satflag, aperture* ap=NULL);

  //x and y are subpixel coordinates
  int wis_photometry(int x, int y, double texp_, int nstack_, vector<double>* phot, int* satflag, aperture* ap=NULL);

  //requires the sub-pixel star position
  //void weighted_photometry(int xx, int yy, double inmag, double* ncounts, double* error, int* satflag, aperture* ap=NULL);

  //requires the sub-pixel star position
  //void ideal_weighted_photometry(int xx, int yy, double inmag, double texp_, int nstack_, double* ncounts, double* error, int* satflag, aperture* ap=NULL);

  //estimate the systematic error due to subpixel shifts from stars within 1
  //pixel of the annulus
  double subpix_systematic(int x, int y, double mag, double frac, starlist* sl, int method);

  //the sky background to the image
  inline void addbg()
  {
    //add background of mag per sq arcsec
    //and any additional background from psf tails
    double Ncounts = bg_counts(1);
    
    for(int i=0;i<Npix;i++)
      {
	timage[i]+=Ncounts;
      }
  };

  //remove the sky background from the image
  inline void subbg()
  {
    //subtract background of mag per sq arcsec
    
    double Ncounts = bg_counts(1);
  
    for(int i=0;i<Npix;i++)
      {
	timage[i]-=Ncounts;
      }
  };

  inline void set_wosrc()
  {
    //copy the current image into memory to act as the knowledge of the 
    //background
  }

  //simulate difference imaging analysis photometry
  //void diaphot(double x, double y, image* ref, double* ncounts, double* error, int* satflag, int ideal=0);
  void diaphot(int x, int y, image* ref, double* incounts, double* ierror, double* ncounts, double* error, int* satflag, int weighted=1);
  void calcrefflux(double x, double y, int weighted=1);
  //create a difference image
  void subimage(image* target, image* ref);

  void setup_fast_photometry(int _xsub, int _ysub, int _x, int _y, double magnitude, double _texp, double _nstack, bool backg=true, aperture* ap=NULL);
  void fast_photometry(double magnification, double* icounts, double* ncounts, double* error, int* satflag, bool ideal=false);

  
 
};
#define IMAGE_HEADER
#endif
