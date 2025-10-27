#include<string>
#include<iostream>
#include<cmath>

#include "psf.h"
#include "image.h"

using namespace std;

// Moffat PSF function (same as in generateMoffat.cpp)
double moffat_psf(double x, double y, void* params)
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
  if(argc < 4 || argc > 5)
    {
      cerr << "Usage: ./generateMoffatPSF <fwhm> <pixel_scale> <output> {<Nkern>}" << endl;
      cerr << "\nOutput file must have a .psf extension" << endl;
      exit(1);
    }

  double fwhm = atof(argv[1]);
  double pixel_scale = atof(argv[2]);
  string output = string(argv[3]);
  int Nkern = 145; // Backwards-compatible default matches historical Roman kernel
  if(argc == 5)
    {
      Nkern = atoi(argv[4]);
      if(Nkern <= 0)
	{
	  cerr << "Nkern must be a positive integer (received " << Nkern << ")" << endl;
	  exit(1);
	}
    }

  // Create PSF with correct subpixel sampling (9x9 grid)
  // Default Nkern=145 matches smoke_test detector configuration (kernelside=291 pix)
  PSF psf(Nkern, 9, pixel_scale);
  
  // Set up PSF parameters
  double beta = 4.0; // Moffat beta parameter
  double params[2] = {fwhm, beta};
  
  cout << "Generating Moffat PSF with subpixel sampling..." << endl;
  cout << "FWHM: " << fwhm << " arcsec, pixel scale: " << pixel_scale << " arcsec/pixel" << endl;
  
  // Generate PSF with subpixel sampling
  // Use the standard step size that the PSF class uses
  double step = pixel_scale / psf.Nsub; // Step size for subpixel sampling
  if(psf.generate_psf(moffat_psf, params, step) < 0)
    {
      cerr << "Failed to generate PSF" << endl;
      exit(1);
    }
  
  // Integrate the PSF
  if(psf.integrate() < 0)
    {
      cerr << "Failed to integrate PSF" << endl;
      exit(1);
    }
  
  cout << "Generated PSF with " << psf.Nsub << "x" << psf.Nsub << " subpixel sampling" << endl;
  
  // Write the binary PSF file
  if(psf.write_psf(output) < 0)
    {
      cerr << "Failed to write PSF file: " << output << endl;
      exit(1);
    }
  
  cout << "PSF written successfully: " << output << endl;
  
  return 0;
}
