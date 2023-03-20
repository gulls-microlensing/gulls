#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<fitsio.h>
#include<ctime>
#include<cstdlib>
#include<string>

#include "psf.h"
#include "image.h"
#include "constants.h"
#include "split.h"
#include "integerPowers.h"

double fs(double Fin, double Fout)
{
  double Fsource = Fin-Fout;
  return Fsource/Fin;
}

int main(int argc, char* argv[])
{
//////////////////////////////////////////////////////////////////////////////
//                               USER OPTIONS
//

  const int maximsize = 1024;            //maximum image size without 
                                         //receiving a warning requiring
                                         //user interaction
  const bool output_image = true;        //output simulated images
  const double large_psf_mag = 15.0;     //magnitude where large psf is used
                                         //missed flux still gets distributed, 
                                         //but uniformly
  const int subpix_method = 0;           //Method by which sub-pixel systematic
                                         //is computed (see below for options)
  const bool zero_systematic = true;     //Do not include a systematic 
                                         //component in the raw error
                                         //can then add either const or 
                                         //calculated systematic

  //Subpixel computation:
  // 0 - the fractional rms fluctuation due the total flux within a slightly 
  //     expanded aperture. (USE THIS)
  // 1 - the fractional rms fluctuation due to the total flux of all stars 
  //     that lie within the same expanded aperture (DO NOT USE THIS)
  //
  //  i.e. method 1 gives a larger error if there are bright stars at the edge 
  //       of the aperture. Method 1 requires significantly more computation
  //       as for each photometered star it must look through all other stars
  //       to see if they are within the aperture. Method 0 just counts up the
  //       total flux in the expanded aperture
                                 
//  
//                               USER OPTIONS
/////////////////////////////////////////////////////////////////////////////

  if(argc!=2)
    {
      cerr << "\nUsage:\n./output_psf <detector>\n" << endl;
      exit(1);
    }

  string detector = string(argv[1]);      //detector parameter file
     
  //Output consists of:
  ////                   <outroot>.txt -the list of stars and photometry
  ////                   <outroot>_true.fits -true starfield image
  ////                   <outroot>_image.fits -single simulated starfield image
  ////                   <outroot>_stack.fits -stacked sim starfield image

  //set up the image

  image im;
  long seed;

  //load the detector - does a lot of the set-up
  if(im.load_detector(detector)<0) exit(1);

  //calculate the systematic error due to sub pixel motions
  //double pixel_systematic = im.sub_pixel_test("subpix.fits");

  int nsub = im.psf.Nsub;
  int nkern = im.psf.Nkern;

  im.set_image_properties(2*nkern+1,2*nkern+1);

  int i=0;

  for(int jsub=0; jsub<1; jsub++)
    {
      for(int isub=0; isub<1; isub++)
	{
	  for(int jk=-nkern;jk<=nkern;jk++)
	    {
	      for(int ik=-nkern;ik<=nkern;ik++)
		{
		  //int i = ik + (2*nkern+1)*(jk + (2*nkern+1)*(isub + nsub*jsub));
		  double val = im.psf.pixval(ik,jk,isub,jsub);
		  im.timage[i] = val;
		  i++;
		}
	    }
	}
    }
  //for(int i=0;i<im.Xpix*im.Ypix;i++) im.timage[i] = im.psf.pixval(i);
  im.write_truefits("PSF_output.fits",true);


  
  im.set_image_properties((2*nkern+1)*nsub,(2*nkern+1)*nsub);

  for(int jsub=0; jsub<nsub; jsub++)
    {
      for(int isub=0; isub<nsub; isub++)
	{
	  for(int jk=-nkern;jk<=nkern;jk++)
	    {
	      for(int ik=-nkern;ik<=nkern;ik++)
		{
		  //int i = ik + (2*nkern+1)*(jk + (2*nkern+1)*(isub + nsub*jsub));
		  int x = (2*nkern+1)*isub + (ik+nkern);
		  int y = (2*nkern+1)*jsub + (jk+nkern);
		  i = x + (2*nkern+1)*nsub*y;
		  double val = im.psf.pixval(ik,jk,isub,jsub);
		  im.timage[i] = val;
		  //i++;
		}
	    }
	}
    }
  //for(int i=0;i<im.Xpix*im.Ypix;i++) im.timage[i] = im.psf.pixval(i);
  im.write_truefits("PSF_output_dith.fits",true);


  for(int jsub=0; jsub<nsub; jsub++)
    {
      for(int isub=0; isub<nsub; isub++)
	{
	  for(int jk=-nkern;jk<=nkern;jk++)
	    {
	      for(int ik=-nkern;ik<=nkern;ik++)
		{
		  //int i = ik + (2*nkern+1)*(jk + (2*nkern+1)*(isub + nsub*jsub));
		  int x = isub + (-ik+nkern)*nsub;
		  int y = jsub + (-jk+nkern)*nsub;
		  i = x + (2*nkern+1)*nsub*y;
		  double val = im.psf.pixval(ik,jk,isub,jsub);
		  im.timage[i] = val;
		  //i++;
		}
	    }
	}
    }
  //for(int i=0;i<im.Xpix*im.Ypix;i++) im.timage[i] = im.psf.pixval(i);
  im.write_truefits("PSF_output_super.fits",true);


  im.set_image_properties(5,5);
  for(i=0;i<5*nsub;i++)
    {
      int x=i;
      int y=int(2.5*nsub);
      im.reset_image();
      im.addstar(x,y,20);
      cout << x << " " << im.timage[11+5] << " " << im.timage[12+5] << " " << im.timage[13+5] << endl;
    }

  

  
}
