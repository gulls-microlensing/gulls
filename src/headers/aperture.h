#ifndef APERTURE_HEADER

#include<vector>
#include<iostream>
#include<cmath>
#include<cstdlib>

using namespace std;

struct aperture
{
  int Naper;
  vector<double> aper;
  bool init;

  aperture()
  {
    init = false;
  }

  void five_nocorners()
  {
    static const int side = 5;  //number of pixels across an appeture - must be odd
    double app[side*side] = {0,1,1,1,0,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  0,1,1,1,0};
    Naper = side;
    aper = vector<double>(app, app+sizeof(app)/sizeof(double));
    init=true;
  }

  void seven_nocorners()
  {
    static const int side = 7;  //number of pixels across an appeture - must be odd
    double app[side*side] = {0,1,1,1,1,1,0,  1,1,1,1,1,1,1,  1,1,1,1,1,1,1,  1,1,1,1,1,1,1,  1,1,1,1,1,1,1,  1,1,1,1,1,1,1,  0,1,1,1,1,1,0};
    Naper = side;
    aper = vector<double>(app, app+sizeof(app)/sizeof(double));
    init=true;
  }

  void square(int side=5)
  {
    Naper = side;
    aper.resize(side*side,1);
    init=true;
  }

  void show_aperture()
  {
    cout << "\nAperture:\n";
    for(int j=0;j<Naper;j++)
      {
	for(int i=0;i<Naper;i++)
	  cout << aper[Naper*j + i] << " ";
	cout << "\n";
      }
    cout << endl;
  }

  //generate a circular aperture
  void generate_aperture(double r, double pixscale)
  {
    r/=pixscale; //work in pixels
    double rad2 = r*r;

    int radius = int(ceil(r-0.5)); //excludes the central pixel
    Naper = 2*radius+1;
    aper.resize(Naper*Naper,1);
    for(int i=0;i<Naper*Naper;i++) aper[i] = 1;

    //work out whether pixel center is within aperture
    double x,y,r2;
    for(int j=0;j<Naper;j++)
      {
	y = (j+0.5) - (radius+0.5);
	for(int i=0;i<Naper;i++)
	  {
	    x = (i+0.5) - (radius+0.5);
	    r2 = x*x + y*y;
	    if(r2>rad2) aper[j*Naper+i] = 0; //outside the circle, set to zero
	  }
      }
    init=true;
  }

  //generate an aperture expanded by one pixel from the aperture that is passed
  void expanded_aperture(aperture* original)
  {
    Naper = original->Naper+2;
    aper.resize(Naper*Naper);

    //work in the original's coordinates
    for(int j=-1;j<=original->Naper;j++)
      {
	for(int i=-1;i<=original->Naper;i++)
	  {
	    //search in the surrounding pixels of the original aperture
	    for(int jj=-1;jj<=1;jj++)
	      {
		int y = j+jj;
		for(int ii=-1;ii<=1;ii++)
		  {
		    int x = i+ii;

		    if(x>=0 && x<original->Naper && y>=0 && y<original->Naper)
		      {
			if(original->aper[j*original->Naper + i])
			  {
			    aper[(j+1)*Naper + (i+1)] = 1;
			  } //if 1
		      } //if in original
		  } //for ii	
	      } //for jj
	  } //for i
      } //for j

    init = true;
  } //function

  //generate a circular aperture
  /*
  void gaussian_aperture(double pixscale, double fwhm)
  {
    double r=pixscale; //work in pixels
    double rad2 = r*r;
    double x,y,r2;
    double sig = 0.5*(fwhm/pixscale)/sqrt(2*log(2));

    int radius = int(ceil(r-0.5)); //excludes the central pixel
    Naper = 2*radius+1;
    aper.resize(Naper*Naper,1);
    for(int j=0;j<Naper;j++) 
      {
	y = (j+0.5) - (radius+0.5);
	for(int i=0;i<Naper;i++)
	  { 
	    x = (i+0.5) - (radius+0.5);
	    r2 = x*x + y*y;
	    //work out whether pixel center is within aperture
	    aper[j*Naper+i] = 1/(2*pi*sig*sig) * exp(-0.5*r2/(sig*sig));
	  }
      }

    for(int j=0;j<Naper;j++)
      {
	y = (j+0.5) - (radius+0.5);
	for(int i=0;i<Naper;i++)
	  {
	    x = (i+0.5) - (radius+0.5);
	    r2 = x*x + y*y;
	    if(r2>rad2) aper[j*Naper+i] = 0; //outside the circle, set to zero
	  }
      }
    init=true;
  }
  */
};

#define APERTURE_HEADER
#endif
