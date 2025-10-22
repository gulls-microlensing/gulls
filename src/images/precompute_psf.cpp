#include<string>
#include<iostream>

#include "psf.h"
#include "image.h"

using namespace std;

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

  if(applyct) 
    {
      im.psf.load_crosstalk(crosstalk);
      im.psf.apply_crosstalk();
    }

  if(im.psf.write_psf(output)); exit(1);

  cout << "PSF written successfully" << endl;
  
}
