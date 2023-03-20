#include<vector>
#include<fstream>
#include<string>
#include<cmath>
#include<iostream>
#include<cstdlib>

#include "split.h"

using namespace std;

int main(int argc, char* argv[])
{

  if(argc!=4)
    {
      cerr << "Usage: ./pixconvolve <input> <kernel> <output>" << endl;
      exit(1);
    }

  string inname = string(argv[1]);
  string kername = string(argv[2]);
  string outname = string(argv[3]);

  ifstream in(inname.c_str());
  if(!in) 
    {
      cerr << "Error: Could not open input " << inname << endl;
      exit(1);
    }
  ifstream ker(kername.c_str());
  if(!ker) 
    {
      cerr << "Error: Could not open kernel " << kername << endl;
      exit(1);
    }
  ofstream out(outname.c_str());
  if(!out) 
    {
      cerr << "Error: Could not open output " << outname << endl;
      exit(1);
    }
  out << "preconvolved" << endl;

  vector<double> input;
  vector<double> kernel;
  int ix=0, iy=0; //input sizes
  int kx=0, ky=0; //kernel sizes
  
  string line;
  vector<double> data;

  //read in the input
  cout << "Reading input..."; cout.flush();
  while(!in.eof())
    {
      getline(in,line);
      split(line,data);

      if(ix==0) 
	{
	  ix = data.size();
	  input.reserve(ix*ix);
	}

      if(int(data.size())==ix)
	{
	  input.insert(input.end(),data.begin(),data.end());
	  iy++;
	}
    }
  cout << "done" << endl;
  cout << "ix iy = " << ix << " " << iy << endl;

  if(ix!=iy)
    {
      cerr << "Warning: mismatched lines in the input file (ix,iy) = (" << ix << "," << iy << ")" << endl;
    }

  //read in the kernel
  while(!ker.eof())
    {
      getline(ker,line);
      split(line,data);

      if(kx==0) 
	{
	  kx = data.size();
	  input.reserve(kx*kx);
	}

      if(int(data.size())==kx)
	{
	  kernel.insert(kernel.end(),data.begin(),data.end());
	  ky++;
	  
	}
    }

  cout << "kx ky = " << kx << " " << ky << endl;

  if(kx!=ky)
    {
      cerr << "Error: mismatched lines in the kernel file (kx,ky) = (" << kx << "," << ky << ")" << endl;
      exit(1);
    }

  //work out the kernel normalization
  double kn=0; //kernel normalization
  for(int i=0;i<int(kernel.size());i++)
    {
      kn+=kernel[i];
    }

  //compute the convolution at the position of each input
  int ii,ji,ik,jk; //input and kernel x and y
  int kshift = (kx-1)/2; //shift the origin of the kernel to the center
  int iidx=0;
  int skip;
  int allskip;
  double sum; //sum

  //for each sampling of the input
  for(ji=0;ji<ix;ji++)
    {
      allskip=1;
      for(ii=0;ii<iy;ii++)
	{
	  ++iidx; //input position index
	  sum=0;
	  skip=0;
	  
	  //for each sampling of the kernel
	  for(jk=0;jk<kx;jk++)
	    {
	      //skip this position if we don't have full input coverage
	      if(ji-kshift<0 || ji+kshift>=iy)
		{
		  skip=1;
		  continue;
		}
		
	      for(ik=0;ik<ky;ik++)
		{
		  //skip this position if we don't have full input coverage
		  if(ii-kshift<0 || ii+kshift>=ix)
		    {
		      skip=1;
		      continue;
		    }

		  sum += input[ix*(ji+jk-kshift)+(ii+ik-kshift)] * kernel[kx*jk+ik];
		}//end for ik
	    } //end for jk

	  if(!skip)
	    {
	      allskip=0;
	      out << sum/kn << " ";
	    }

	} //end for ii
      if(!allskip) out << endl;
      cout << "."; cout.flush();
    } //end for ji
  cout << endl;	  
}

