#include<string>
#include<iostream>

#include "psf.h"
#include "image.h"

using namespace std;

// Moffat PSF function for generating PSFs from detector parameters
double eval_moffat(double x, double y, void* params)
{
  double* p = (double*)params;
  double fwhm = p[0];
  double beta = p[1];
  
  double r = sqrt(x*x + y*y);
  double alpha = fwhm / (2.0 * sqrt(pow(2.0, 1.0/beta) - 1.0));
  
  return pow(1.0 + (r/alpha)*(r/alpha), -beta);
}

int main(int argc, char* argv[])
{
  if(argc!=3&&argc!=4)
    {
      cerr << "Usage: ./precompute_psf <detector> {<crosstalk>} <output>" << endl;
      cerr << "\nOutput file must have a .psf extension" << endl;
      exit(1);
    }

  string detector = string(argv[1]);
  string crosstalk; int applyct=0;
  string output = string(argv[argc-1]);

  if(argc==4) 
    {
      applyct=1;
      crosstalk = string(argv[2]);
    }

  image im;

  if(im.load_detector(detector)<0) exit(1);
  cout << "detector read" << endl;

  // Generate PSF from detector parameters if no PSF file was loaded
  if(!im.psf.init)
    {
      cout << "Generating PSF from detector parameters..." << endl;
      
      // Get PSF parameters from detector
      double fwhm = im.psf.pixscale * 0.2; // PSFFWHM in arcsec, convert to pixel scale
      double beta = 4.0; // Moffat beta parameter
      
      double params[2] = {fwhm, beta};
      
      // Generate PSF with subpixel sampling
      double step = im.psf.pixscale / im.psf.Nsub; // Step size for subpixel sampling
      if(im.psf.generate_psf(eval_moffat, params, step) < 0)
        {
          cerr << "PSF generation failed" << endl;
          exit(1);
        }
      
      if(im.psf.integrate() < 0)
        {
          cerr << "PSF integration failed" << endl;
          exit(1);
        }
      
      cout << "PSF generated successfully" << endl;
    }

  if(applyct) 
    {
      im.psf.load_crosstalk(crosstalk);
      im.psf.apply_crosstalk();
    }

  if(im.psf.write_psf(output)) exit(1);

  cout << "PSF written successfully" << endl;
  
}
