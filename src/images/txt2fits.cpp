#include<string>
#include<iostream>

#include "psf.h"
#include "image.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=3&&argc!=4)
    {
      cerr << "Usage: ./precompute_psf filename <output>" << endl;
      cerr << "\nOutput file must have a .fits extension" << endl;
      exit(1);
    }

  string input = string(argv[1]);
  string output = string(argv[argc-1]);

  PSF inpsf;

  inpsf.load_txt(input, 1.0/11.0);

  


  if(inpsf.write_psf(output)) exit(1);


  
}
