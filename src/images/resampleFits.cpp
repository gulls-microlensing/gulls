#include<string>
#include<iostream>

#include<fitsio.h>

using namespace std;

int main(int argc, char* argv[])
{

  if(argc!=4&&argc!=5)
    {
      cerr << "Usage: ./resampleFits <input> <resampling> {<flux scale>} <output>" << endl;
      exit(1);
    }

  fitsfile *out;
  fitsfile *in;

  string infile = string(argv[1]);
  int fac = atoi(argv[2]); //resolution multiplication factor
  string outfile = string(argv[argc-1]);

  int inaxis=2;
  long inaxes[2] = {1,1};
  long onaxes[2] = {1,1};
  long firstpix[2] = {1,1};
  long ofirstpix[2] = {1,1};

  double fluxscale=1;
  if(argc==5) fluxscale=atof(argv[3]);

  int status=0;

  double* idata;
  double* odata;

  fits_open_file(&in, infile.c_str(), READONLY, &status);

  fits_get_img_dim(in, &inaxis, &status);
  fits_get_img_size(in, inaxis, inaxes, &status);

  onaxes[0] = abs(fac)*inaxes[0];
  onaxes[1] = abs(fac)*inaxes[1];

  idata = new double[inaxes[0]];
  odata = new double[onaxes[0]];

  if(fac>0)
    {
      fits_create_file(&out,outfile.c_str(),&status);
      fits_create_img(out, DOUBLE_IMG, inaxis, onaxes, &status);
      
      int nrow = inaxes[0]; //number of pixels in a row

      //loop over the rows
      for(firstpix[1]=1; firstpix[1]<=inaxes[1]; firstpix[1]++)
	{
	  fits_read_pix(in, TDOUBLE, firstpix, nrow, NULL, idata, NULL, 
			&status);

	  //loop over the output rows
	  for(ofirstpix[1] = (firstpix[1]-1)*fac+1; ofirstpix[1] <= (firstpix[1]-1)*fac+3; ofirstpix[1]++)
	    {
	      for(int i=0;i<nrow;i++)
		{
		  for(int j=0;j<fac;j++)
		    {
		      odata[i*fac+j] = fluxscale*idata[i];
		    }
		}
	      
	      fits_write_pix(out, TDOUBLE, ofirstpix, fac*nrow, odata, 
			     &status);

	    }	  
	}
    }
  else
    {
      cerr << "Error: reducing the resolution is not currently supported\n";
      exit(1);
    }

  fits_close_file(in,&status);
  fits_close_file(out,&status);

  delete[] idata;
  delete[] odata;

}
