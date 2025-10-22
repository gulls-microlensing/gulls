#include<string>
#include<iostream>

#include "psf.h"
#include "image.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=4)
    {
      cerr << "Usage: ./generate_psf_from_text <text_psf> <pixel_scale> <output>" << endl;
      cerr << "\nOutput file must have a .psf extension" << endl;
      exit(1);
    }

  string input = string(argv[1]);
  double pixel_scale = atof(argv[2]);
  string output = string(argv[3]);

  PSF psf;
  
  // Load the text PSF file
  if(psf.load_txt(input, pixel_scale) < 0)
    {
      cerr << "Failed to load text PSF file: " << input << endl;
      exit(1);
    }
  
  cout << "Loaded text PSF with " << psf.Nsub << "x" << psf.Nsub << " subpixel sampling" << endl;
  
  // Write the binary PSF file
  if(psf.write_psf(output) < 0)
    {
      cerr << "Failed to write PSF file: " << output << endl;
      exit(1);
    }
  
  cout << "PSF written successfully: " << output << endl;
  
  return 0;
}
