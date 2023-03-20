#include<vector>
#include<cmath>
#include<string>
#include<fstream>
#include<iostream>
#include<fitsio.h>
#include<cstdlib>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>

#include "psf.h"
#include "constants.h"
#include "integerPowers.h"
#include "split.h"

using namespace std;

/**********************************************************************

                               image + psf

                        Written by: Matthew Penny
                         Department of Astronomy
                          Ohio State University
                               Formerly
                 Jodrell Bank Centre for Astrophysics
                        University of Manchester

                      Copyright: Matthew Penny 2011-2012

A C++ class for simulating imaging and photometry in dense star-fields

If you wish to use this software for a publication, please contact me
before doing too much work, as I may require to be an author.
                     

**********************************************************************/

double eval_gaussian(double x, double y, void* params)
{
  double fwhm = *(double *) params;
  double s2=sqr(fwhm)/(8*log(2));
  double r2=x*x+y*y;
  double val=exp(-0.5*r2/s2)/(twoPi*s2); //2D normalized gaussian 
  return val;
}
/*
void PSF::integrate(double (*psffunc)(double, double, void*), void* params)
{
  //Outputs a table of the PSF integrated over a pixel at sub pixel intervals

  //Nkern is the half kernel size

  vector<double> points = vector<double>(sqr(pointsside));
  psf.resize(psfsize);

  int i, j;
  double x,y;

  //calculate the psf value at discrete points
  for(j=0; j<pointsside; j++)
    {
      y = double(j-origin)*pixscale/double(Nsub);
      for(i=0; i<pointsside; i++)
	{
	  x = double(i-origin)*pixscale/double(Nsub);
	  //if(j==0&&i<Nsub) cout << i << " " << x << endl;
	  points[pointsside*j+i] = psffunc(x,y,params);
	}
    }
  //cout << "Points calculated" << endl;

  //integrate the psf over pixels

  //if a star is positioned at xs,ys within a pixel,
  //its pixel integrated psf will be an Nkern*Nkern array starting at record
  // sqr(1+2*(Nkern-1))*(Nsub*ys+xs) - i.e. there is a kernel for each sub
  //pixel star placement stored in contiguous array memory. The kernel runs
  //over the pixels from -(Nkern-1) to (Nkern+1) surrounding the pixel
  //containing the star center

  const double ot=1.0/3.0;

  //for each y-sub-pixel star location
  for(int ys=0;ys<Nsub;ys++)
    {
      //for each x-sub-pixel point
      for(int xs=0;xs<Nsub;xs++)
	{
	  //the position of the first pixel of the kernel in the array
	  int astart = kernsize*(Nsub*ys+xs);

	  //for each y-pixel in the kernel
	  for(int yp=-Nkern;yp<=Nkern;yp++)
	    {
	      
	      //the position of the integration limits are:
	      //         sub pixel offset   kernel pixel offset
	      int ymin = -ys              + (yp*Nsub) + origin;
	      int ymax = -ys + (Nsub)     + (yp*Nsub) + origin;

	      //	      if(xs==0 && ys==0) cout << "y: " << yp << " " << ymin << " " << ymax << " " << double(ymin-origin)*pixscale/double(Nsub) << " " << double(ymax-origin)*pixscale/double(Nsub) <<  endl;

	      //for each x-pixel in the kernel
	      for(int xp=-Nkern;xp<=Nkern;xp++)
		{
		  //the position of the pixel in the array
		  int ashift = kernside*(yp+Nkern) + (xp+Nkern);
		  int apos = astart+ashift;

		  //integrate the psf over this pixel

		  //the position of the integration limits are:
		  //         sub pixel offset   kernel pixel offset
		  int xmin = -xs              + (xp*Nsub) + origin;
		  int xmax = -xs + (Nsub)     + (xp*Nsub) + origin;

		  //		  if(xs==0 && ys==0) cout << "x: " << xp << " " << xmin << " " << xmax << " " << double(xmin-origin)*pixscale/double(Nsub) << " " << double(xmax-origin)*pixscale/double(Nsub) << " " << points[pointsside*ymin+xmin] << " " << points[pointsside*ymin+xmax] << endl;

		  //integrate over the limits
		  psf[apos] = xint(xmin,xmax,ymin,&points);
		  for(int j=ymin+1; j<ymax-(Nsub%2?1:0); j++)
		    {
		      psf[apos] += ((j-ymin)%2?4.0:2.0) 
			* xint(xmin,xmax,j,&points);
		    }
		  double ep = xint(xmin,xmax,ymax-(Nsub%2?1:0),&points);
		  if(Nsub%2) //if we have an odd no of points, tack on a trapezium
		    psf[apos] += 2.5*ep + 1.5*xint(xmin,xmax,ymax,&points);
		  else psf[apos]+=ep;
		  
		  psf[apos] *= sqr(ot*pixscale/double(Nsub));
		      
		} //end for each x-pixel in the kernel

	    } //end for each y-pixel in the kernel

	} //end for each x-sub-pixel point

    } //end for each y-sub-pixel star location

  //The psf is now safe to use
  initialized=true;
  //  cout << "Integrated" << endl;

}
*/

/*
double PSF::xint(int xmin, int xmax, int y, vector<double>* p)
{
  //integrate the psf along the x-axis using simpson's extended rule
  //use abs to get values in other quadrants

  //cout << "enter xint" << endl;

  int offset = pointsside*y;

  //cout << offset + xmin << " " << offset+xmax << " " << p->size() << endl;

  double sum = (*p)[offset+xmin];
  for(int i=xmin+1; i<xmax-(Nsub%2?1:0); i++)
    {
      sum += ((i-xmin)%2?4.0:2.0) * (*p)[offset+i];
    }
  double ep = (*p)[offset+xmax-(Nsub%2?1:0)];
  if(Nsub%2) //if we have an odd no of points, tack on a trapezium
    sum += 2.5*ep + 1.5*(*p)[offset+xmax];
  else sum+=ep;

  //cout << "exit xint" << endl;

  return sum;
}
*/

/*
double PSF::xint(int xmin, int xmax, int y, double p[])
{
  //integrate the psf along the x-axis using simpson's extended rule
  //use abs to get values in other quadrants

  //cout << "enter xint" << endl;

  int offset = pointsside*y;

  //cout << offset + xmin << " " << offset+xmax << " " << p->size() << endl;

  double sum = p[offset+xmin];
  for(int i=xmin+1; i<xmax-(Nsub%2?1:0); i++)
    {
      sum += ((i-xmin)%2?4.0:2.0) * p[offset+i];
    }
  double ep = p[offset+xmax-(Nsub%2?1:0)];
  if(Nsub%2) //if we have an odd no of points, tack on a trapezium
    sum += 2.5*ep + 1.5*p[offset+xmax];
  else sum+=ep;

  //cout << "exit xint" << endl;

  return sum;
}*/

void PSF::output_psf(int si, int sj, string filename, bool overwrite)
{
 
  fitsfile *out;

  int status=0;
  long fpixel=1;
  long naxis=2;
  long naxes[2] = {kernside,kernside};
  long nelements=kernside*kernside;

  //delete the old file if overwrite is true
  if(overwrite)
    {
      ifstream tempfile(filename.c_str());
      if(tempfile)
	{
	  tempfile.close();
	  if(remove(filename.c_str())) 
	    {
	      perror("Error deleting file - will not try writing fits as it will probably fail anyway.");
	      return;
	    }
	}
    }

  //temporarily grab some memory
  float *array = new float[nelements];

  //write the array
  for(int j=0;j<kernside;j++)
    {
      for(int i=0;i<kernside;i++)
	{      
	  array[kernside*j+i] = psf[(Nsub*sj+si)*kernsize + kernside*j+i];
	}
    }

  fits_create_file(&out,filename.c_str(),&status);
  fits_create_img(out, FLOAT_IMG, naxis, naxes, &status);

  fits_write_img(out, TFLOAT, fpixel, nelements, array, &status);

  fits_close_file(out, &status);
  fits_report_error(stderr, status);

  //remember to delete the memory
  delete [] array;

  //outputs in a different convention to that in which its stored

  cout << "Nsub = " << Nsub << endl;
  cout << "Nkern = " << Nkern << endl;
  cout << "psf.size() = " << psf.size() << endl;
  cout << "psf.size(), /kernsize,  = " << psf.size() << " " << psf.size()/kernsize << endl;

  /*  for(int j=-Nkern;j<=-Nkern;j++)
    {
      for(int i=-Nkern;i<=Nkern;i++)
	{
	  out << i << " " << j << " " << pixval(si,sj,i,j) << "\n";
	}
      out << "\n";
      }*/
}

//Load a psf from a text file
int PSF::load_txt(string fname, double step, double centx, double centy, int nx, int ny)
{
  vector<double> storage;
  vector<double> data;
  string line;
  int npix;

  psfh = step;
  bool autocentre=false;

  if(centx==1e30 || centy==1e30) autocentre=true;

  double maxpix=-1e99;
  int centpix=-1;

  ifstream in(fname.c_str());
  if(!in)
    {
      cerr << __FUNCTION__ << "Error: Could not open text file (" << fname << ")" << endl;
      exit(1);
    }

  int firstline=1;
  preconvolved=0;

  while(!in.eof())
    {
      getline(in,line);

      if(firstline)
	{
	  if(line.find("preconvolved")!=string::npos)
	    {
	      cout << "PSF has been preconvolved - will not integrate" << endl;
	      preconvolved=1;
	      continue;
	    }
	}
      split(line,data);

      for(int i=0;i<int(data.size());i++)
	storage.push_back(data[i]);
    }

  npix = nx*ny;

  if(nx<=0 ||  ny<=0)
    {
      //work out how big the image is - assume square
      npix = storage.size();

      if(npix % int(sqrt(npix)) !=0 && int(npix / sqrt(npix)) != npix)
	{
	  //image isn't square
	  cerr << __FUNCTION__ << ": Error: Tried to load a non-square image - if this was intended, specify the number of pixels on each side." << endl;
	  cerr << "npix = " << npix << endl;
	  exit(1);
	}

      nx = ny = int(sqrt(npix));
    }

  if(int(storage.size()) != npix)
    {
      cerr << __FUNCTION__ << ": Error: Number of pixels is not what is expected:" << endl;
      cerr << "(nx, ny, npix, nactualpix) = (" << nx << ", " << ny << ", " << npix << ", " << storage.size() << endl;
      exit(1);
    }

  psfx = nx;
  psfy = ny;

  //check if there is already an array loaded
  if(psfxy!=NULL)
    {
      delete[] psfxy;
      psfxy=NULL;
    }

  //reserve space
  psfxy = new double[psfx*psfy];

  for(int i=0;i<psfx*psfy;i++)
    {
      psfxy[i] = storage[i];
      if(autocentre)
	{
	  if(psfxy[i]>maxpix)
	    {
	      maxpix = psfxy[i];
	      centpix = i;
	    }
	}
    }

  //work out the actual coordinate of the bottom left corner of the grid
  if(autocentre)
    {
      if(centpix<0 || centpix>=psfx*psfy)
	{
	  cerr << __FUNCTION__ << ": Warning: A peak pixel was not found - psf will not be centered" << endl;
	}
      ox = -(centpix%psfx)*psfh;
      oy = -floor(centpix/psfx)*psfh;
    }
  else
    {
      ox = -0.5*psfx*psfh;
      oy = -0.5*psfy*psfh;
    }

  //verification
  cout << "Verification:" << endl;
  cout << "pixsize (\"): " << pixscale << endl;
  cout << "(ox,oy) (\"): " << "(" << ox << "," << oy << ")" << endl;
  cout << "(psfx, psfy) = (" << psfx << ", " << psfy << ")" << endl;
  cout << "PSF centre(\"): x = " << -0.5*psfx*psfh - ox << ", y = " 
       << -0.5*psfy*psfh - oy << endl;
  cout << "PSF centre(txtfile pixels): x = " << -0.5*psfx - ox/psfh << ", y = " 
       << -0.5*psfy - oy/psfh << endl;
  cout << "Corner vals:" << endl;
  cout << psfxy[nx*(ny-1) + 0] << "\t" << psfxy[nx*(ny-1) + (nx-1)] << endl;
  cout << psfxy[nx*0 + 0] << "\t" << psfxy[nx*0 + (nx-1)] << endl;
  cout << endl;

  //interpolate
  interpolate();

  return 0;
}

//generate a psf given a function
int PSF::generate_psf(double (*psffunc)(double, double, void*), void* params, double step)
{
  double x,y;

  init = false;

  //the step size
  psfh = step;

  //we can be efficient and calculate the psf only over those points we need
  //needs to cover at least +1 kernel pixel either side

  //number of samples to a side
  psfx = psfy = long((kernside+2)*ceil(pixscale/psfh));
  //cerr << "ceil(pixscale/psfh) = " << ceil(pixscale/psfh) << endl;

  //allocate memory
  psfxy = new double[psfx*psfy];

  //poisition of the first point in coords centered on the psf center
  ox = -0.5*psfx*psfh; 
  oy = -0.5*psfy*psfh; 

  //calculate the psf value at discrete points
  for(int j=0; j<psfy; j++)
    {
      y = oy + double(j)*psfh;
      for(int i=0; i<psfx; i++)
	{
	  x = ox + double(i)*psfh;
	  psfxy[psfx*j+i] = psffunc(x,y,params);
	}
    }

  //compute a set of spline interpolations
  interpolate();

  return 0;
}

/*
double PSF::integrate_pixel(double x0, double y0, double x1, double y1)
{
  //integrate the tabulated psf over the pixel defined by these corners
  //origin for the input is bottom left of central pixel input is in units of
  //arcsec

  double oy = -0.5*psfy*psfh; //position of the input origin in coordinates of 
                              //the tabulated psf

  int f = int(ceil((y0 - oy)/psfh)); //index of first integration point
  double fy = f*psfh + oy; //actual position of the first point
  int l = int(floor((y1 - oy)/psfh)); //index of last integration point
  double ly = l*psfh + oy; //actual position of the first point

  double sum;

  double r0 = integrate_x(x0,x1,f);
  double r1 = integrate_x(x0,x1,l);

  sum = 0.5*(r0 + r1);
  for(int j=f+1;j<l;j++)
    {
      sum+= integrate_x(x0,x1,f);
    }
  sum*=psfh;

  //+ funny endpoints to cope with arbitrary placement
  //trap rule over a small inteval using the interpolated intensity
  //sum += 0.5 * (fy-y0) 
  //  * (r0+(y0-(fy-psfh))*(r0+integrate_x(x0,x1,f-1))/psfh);
  //sum += 0.5 * (y1-ly) 
  //  * (r1+(y1-ly)*(integrate_x(x0,x1,l+1)+r1)/psfh);

  return sum;

}
*/

/*
double PSF::integrate_x(double x0, double x1, int y)
{  
  //y is the row number in the table of points
  int shift = y*psfx;

  double ox = -0.5*psfx*psfh; //position of point i=0

  int f = int(ceil((x0 - ox)/psfh)); //index of first integration point
  double fx = f*psfh + ox; //actual position of the first point
  int l = int(floor((x1 - ox)/psfh)); //index of last integration point
  double lx = l*psfh + ox; //actual position of the last point

  double sum;
  //cout << "(x0,x1,y,ox,f,l) = (" << x0 << ", " << x1 << ", " << y << ", " << ox << ", " << f << ", " << l << ")" << endl; 
  //cout << "(psfx, psfy, shift, f, l, psfx*psfy) = (" << psfx << ", " << psfy << ", " << shift << ", " << f << ", " << l << ", " << psfx*psfy << ")" << endl; 
  //trapezoidal rule
  sum = 0.5*(psfxy[shift+f] + psfxy[shift+l]);
  //if(x0==0) cout << "(fx,lx,y) = (" << fx << ", " << lx << ", " << y*psfh-0.5*psfy*psfh << ")" << endl;
  //if(x0==0) cout << psfxy[shift+f] << " ";
  for(int i=f+1; i<l; i++)
    {
      sum += psfxy[shift+i];
      //if(x0==0.0) cout << psfxy[shift+i] << " ";
    }
  //if(x0==0.0) cout << psfxy[shift+l] << endl;
  sum*=psfh;

  //+ funny endpoints to cope with arbitrary placement
  //trap rule over a small inteval using the interpolated intensity
  //sum += 0.5 * (fx-x0) 
  //  * (psfxy[shift+f]+(x0-(fx-psfh))*(psfxy[shift+f]+psfxy[shift+f-1])/psfh);
  //sum += 0.5 * (x1-lx) 
  //  * (psfxy[shift+l]+(x1-lx)*(psfxy[shift+l+1]+psfxy[shift+l])/psfh);

  //if(x0==0)
  //cout << sum << " " << sum - 0.5 * (x1-lx) 
  //  * (psfxy[shift+l]+(x1-lx)*(psfxy[shift+l+1]+psfxy[shift+l])/psfh) - 0.5 * (fx-x0) 
  //  * (psfxy[shift+f]+(x0-(fx-psfh))*(psfxy[shift+f]+psfxy[shift+f-1])/psfh) << " " << 0.5 * (fx-x0) 
  //  * (psfxy[shift+f]+(x0-(fx-psfh))*(psfxy[shift+f]+psfxy[shift+f-1])/psfh) << " " << 0.5 * (x1-lx) 
  //  * (psfxy[shift+l]+(x1-lx)*(psfxy[shift+l+1]+psfxy[shift+l])/psfh) << endl;
    
    return sum;
}
*/

//integrate a psf that has been calculated/loaded
int PSF::integrate()
{
  //check to see if the psf is available
  if(psfxy==NULL)
    {
      cerr << __FUNCTION__ << ": Error: psf has not been calculated or loaded" << endl;
      return 1;
    }

  cerr << "psfsize = " << psfsize << endl;

  //allocate space
  psf.resize(psfsize);
  spsf.resize(spsfsize);

  //for each y-sub-pixel star location
  for(int ys=0;ys<Nsub;ys++)
    {
      //for each x-sub-pixel point
      for(int xs=0;xs<Nsub;xs++)
	{
	  if(!preconvolved) integrate_kernel(-1,-1,xs,ys);
	  else extract_kernel(-1,-1,xs,ys);
	} //end for each x-sub-pixel point

    } //end for each y-sub-pixel star location

  init = true;
  cout << endl;

  return 0;

}

//integrate the psf centered at (x,y) over the pixels of a kernel 
void PSF::integrate_kernel(double x, double y, int xs, int ys)
{
  //xs,ys are the sub-pixel coordinates of the kernel center
  //set to -1,-1 if not calculating grid of kernels
  //if xs,ys are valid, x,y is ignored
  //x,y are in units of pixels

  int astart;
  int sstart;
  vector<double>* target;

  double sum = 0;
  double ssum = 0;

  if(xs<-1 || xs>=Nsub || ys<-1 || ys>=Nsub)
    {
      cerr << __FUNCTION__ << ": Error: invalid subpixel coordinates (" << xs << "," << ys << ")" << endl;
      exit(1);
    }

  if(xs==-1||ys==-1)
    {
      astart=0;
      sstart=0;
      ipsf.resize(kernsize);
      target = &ipsf;
    }
  else
    {
      x = double(xs)/double(Nsub);
      y = double(ys)/double(Nsub);
      //the position of the first pixel of the kernel in the array
      astart = kernsize*(Nsub*ys+xs);
      sstart = skernsize*(Nsub*ys+xs);
      target = &psf;
    }

  printf("Integrating kernel at (x,y) = (%f,%f)\r",x,y);
  fflush(stdout);

  //cout << "sstart = " << sstart << endl;

  //check that there is enough psf data to cover the kernel

  //needs generalizing to general psf center
  //  if((-Nkern-Nsub)*pixscale < max(ox,oy) || (Nkern+Nsub)*pixscale > min(ox+psfx*psfh,oy+psfy*psfh)) - this is almost certainly wrong...
  if((-Nkern-1)*pixscale < max(ox,oy) || (Nkern+1)*pixscale > min(ox+psfx*psfh,oy+psfy*psfh)) //the one is to allow subpixel movement
    {
      cerr << __FUNCTION__ << ": Error: There is not enough data in the psf that is loaded:" << endl;
      //cerr << "Kernel covers " << (-Nkern-Nsub)*pixscale << " < x < " << (Nkern+Nsub)*pixscale << "; " << (-Nkern-Nsub)*pixscale << " < y < " << (Nkern+Nsub)*pixscale << endl;
      cerr << "Kernel covers " << (-Nkern-1)*pixscale << " < x < " << (Nkern+1)*pixscale << "; " << (-Nkern-1)*pixscale << " < y < " << (Nkern+1)*pixscale << endl;
      cerr << "Tabulated PSF covers " << ox << " < x < " << ox+psfx*psfh << "; " << oy << " < y < " << oy+psfy*psfh << endl;
      exit(1);
    }
  
  //for each y-pixel in the kernel
  for(int yp=-Nkern;yp<=Nkern;yp++)
    {
      //the position of the integration limits are:
      double ymin = (yp - y) * pixscale;
      double ymax = ymin + pixscale;
      
      //for each x-pixel in the kernel
      for(int xp=-Nkern;xp<=Nkern;xp++)
	{
	  //the position of the integration limits are:
	  double xmin = (xp - x) * pixscale;
	  double xmax = xmin + pixscale;

	  //if(ys==0&& xs==0) cout << "ys = 0: " << xmin << " " << xmax << endl;

	  //the position of the pixel in the array
	  int ashift = kernside*(yp+Nkern) + (xp+Nkern);
	  int apos = astart+ashift;

	  double pixval = integrate_spline(xmin,ymin,xmax,ymax);

	  //if(xs==0 && ys==0)
	  //  cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << pixval << endl;

	  (*target)[apos] = pixval;
	  sum+=pixval;

	  if(xs!=-1 || ys!=-1)
	    {
	      if(xp>=-sNkern && xp<=sNkern && yp>=-sNkern && yp<=sNkern)
		{
		  //add to the smaller kernel as well
		  //cout << "sidx = " << sstart + skernside*(yp+sNkern) + (xp+sNkern) << " " << spsf.size() << endl;
		  spsf[sstart + skernside*(yp+sNkern) + (xp+sNkern)] += pixval;
		  //ssum+=pixval;
		}
	    }
	} //end for each x-pixel
      //if(xs==0&&ys==0) cout << endl;
      
    } //end for each y-pixel
  
  //normalize the kernel - assumes all flux within kernel
  for(int i=astart;i<astart+kernsize;i++)
    {
      (*target)[i] /= sum;
    }

  for(int i=sstart;i<sstart+skernsize;i++)
    {
      spsf[i] /= sum;
      ssum += spsf[i];
    }
  //keep track of the missing flux
  //cout << sum << " " << ssum << endl;
  sMissingFlux += (1 - ssum)/sqr(Nsub);
}

//Build interpolation functions for the PSF
void PSF::interpolate()
{
  double* x = new double[psfx];
  double* y = new double[psfy];

  //for each row
  for(int j=0;j<psfy;j++)
    {
      int shift = j*psfx;
      //assign values to the arrays
      for(int i=0;i<psfx;i++)
	{
	  x[i] = ox + double(i)*psfh;
	  y[i] = psfxy[shift + i];
	}

      //cout << "spline size " << psfx << endl;
      acc.push_back(gsl_interp_accel_alloc());
      spline.push_back(gsl_spline_alloc(gsl_interp_cspline,psfx));
      //cout << "spline allocated" << endl;

      gsl_spline_init(spline[j],x,y,psfx);
      //cout << "spline calculated" << endl;

    }

  delete[] x;
  delete[] y;
}

double PSF::integrate_spline(double x0, double y0, double x1, double y1)
{
  //integrate the spline interpolation of the psf over a pixel

  //will need to interpolate over the rows as well, so set up another 
  //interpolation

  //first work out the row numbers of the rows we need to integrate
  int r0 = int(floor((y0 - oy)/psfh));
  int r1 = int(ceil((y1 - oy)/psfh));
  int nr = r1-r0+1;

  double* x = new double[nr]; //actually y values
  double* y = new double[nr]; //actually int_x0^x1 psfxy dx

  //cout << "spline: from " << r0 << " to " << r1 << ", size = " << r1-r0+1 << endl;

  gsl_interp_accel *racc = gsl_interp_accel_alloc();
  gsl_spline *rspline = gsl_spline_alloc(gsl_interp_cspline,nr);

  //calculate the x and y values
  for(int j=r0;j<=r1;j++)
    {
      //cout << "j, j-r0 = " << j << ", " << j-r0 << endl;
      x[j-r0] = oy + j*psfh;
      y[j-r0] = gsl_spline_eval_integ(spline[j], x0, x1, acc[j]);
      //y[j-r0] = 1;
    }

  //now create the spline and integrate it
  gsl_spline_init(rspline, x, y, r1-r0+1);
  double result = gsl_spline_eval_integ(rspline, y0, y1, racc);

  //free up memory
  gsl_spline_free(rspline);
  gsl_interp_accel_free(racc);

  flush_psfxy(); //we have no need for the tabulated psf now

  delete[] x;
  delete[] y;

  //return the result
  return result;
}

//extract the psf value at the center (x,y) of the pixels of a kernel 
void PSF::extract_kernel(double x, double y, int xs, int ys)
{
  //xs,ys are the sub-pixel coordinates of the kernel center
  //set to -1,-1 if not calculating grid of kernels
  //if xs,ys are valid, x,y is ignored
  //x,y are in units of pixels

  int astart;
  int sstart;
  vector<double>* target;

  double sum = 0;
  double ssum = 0;

  if(xs<-1 || xs>=Nsub || ys<-1 || ys>=Nsub)
    {
      cerr << __FUNCTION__ << ": Error: invalid subpixel coordinates (" << xs << "," << ys << ")" << endl;
      exit(1);
    }


  if(xs==-1||ys==-1)
    {
      astart=0;
      sstart=0;
      ipsf.resize(kernsize);
      target = &ipsf;
    }
  else
    {
      x = double(xs)/double(Nsub);
      y = double(ys)/double(Nsub);
      //the position of the first pixel of the kernel in the array
      astart = kernsize*(Nsub*ys+xs);
      sstart = skernsize*(Nsub*ys+xs);
      target = &psf;
    }

  //cout << "sstart = " << sstart << endl;

  //check that there is enough psf data to cover the kernel

  //needs generalizing to general psf center
  //  if((-Nkern-Nsub)*pixscale < max(ox,oy) || (Nkern+Nsub)*pixscale > min(ox+psfx*psfh,oy+psfy*psfh))
  if((-Nkern-1)*pixscale < max(ox,oy) || (Nkern+1)*pixscale > min(ox+psfx*psfh,oy+psfy*psfh))
    {
      cerr << __FUNCTION__ << ": Error: There is not enough data in the psf that is loaded:" << endl;
      cerr << "Kernel covers " << (-Nkern-1)*pixscale << " < x < " << (Nkern+1)*pixscale << "; " << (-Nkern-1)*pixscale << " < y < " << (Nkern+1)*pixscale << endl;
      cerr << "Tabulated PSF covers " << ox << " < x < " << ox+psfx*psfh << "; " << oy << " < y < " << oy+psfy*psfh << endl;
      exit(1);
    }
  
  //for each y-pixel in the kernel
  for(int yp=-Nkern;yp<=Nkern;yp++)
    {
      //the position of the integration limits are:
      double ymin = (yp - y + 0.5) * pixscale;
      
      //for each x-pixel in the kernel
      for(int xp=-Nkern;xp<=Nkern;xp++)
	{
	  //the position of the integration limits are:
	  double xmin = (xp - x + 0.5) * pixscale;

	  //the position of the pixel in the array
	  int ashift = kernside*(yp+Nkern) + (xp+Nkern);
	  int apos = astart+ashift;

	  double pixval = extract_spline(xmin,ymin);

	  (*target)[apos] = pixval;
	  sum+=pixval;

	  if(xs!=-1 || ys!=-1)
	    {
	      if(xp>=-sNkern && xp<=sNkern && yp>=-sNkern && yp<=sNkern)
		{
		  //add to the smaller kernel as well
		  spsf[sstart + skernside*(yp+sNkern) + (xp+sNkern)] += pixval;
		}
	    }
	} //end for each x-pixel
      
    } //end for each y-pixel
  
  //normalize the kernel - assumes all flux within kernel
  for(int i=astart;i<astart+kernsize;i++)
    {
      (*target)[i] /= sum;
    }

  for(int i=sstart;i<sstart+skernsize;i++)
    {
      spsf[i] /= sum;
      ssum += spsf[i];
    }
  //keep track of the missing flux
  sMissingFlux += (1 - ssum)/sqr(Nsub);
}

double PSF::extract_spline(double x0, double y0)
{
  //extract the spline interpolation of the psf at the centre of a pixel

  //will need to interpolate over the rows as well, so set up another 
  //interpolation

  //first work out the row numbers of the rows we need to integrate
  int r0 = int(floor((y0 - 2.5*psfh - oy)/psfh));
  int r1 = int(ceil((y0 + 2.5*psfh - oy)/psfh));
  int nr = r1-r0+1;

  double* x = new double[nr]; //actually y values
  double* y = new double[nr]; //actually psfxy

  gsl_interp_accel *racc = gsl_interp_accel_alloc();
  gsl_spline *rspline = gsl_spline_alloc(gsl_interp_cspline,nr);

  //calculate the x and y values
  for(int j=r0;j<=r1;j++)
    {
      x[j-r0] = oy + j*psfh;
      y[j-r0] = gsl_spline_eval(spline[j], x0, acc[j]);
    }

  //now create the spline and integrate it
  gsl_spline_init(rspline, x, y, nr);
  double result = gsl_spline_eval(rspline, y0, racc);

  //free up memory
  gsl_spline_free(rspline);
  gsl_interp_accel_free(racc);

  flush_psfxy(); //we have no need for the tabulated psf now

  delete[] x;
  delete[] y;

  //return the result
  return result;
}

//write out the psfs to a binary file
int PSF::write_psf(string filename)
{
  int status = 0;

  ofstream out(filename.c_str(),ios::out | ios_base::binary);
  if(!out)
    {
      cerr << __FUNCTION__ << ": Error: Could not open splines file for output (" << filename << ")" << endl;
      status = 1;
      return status;
    }

  cout << "file open" << endl;

  //write out the information needed to use the splines

  //need to be hard coded by the developer. Any time a change is made increment
  //the version num counter (by hand, here, and in the write function) so that
  //file format conflicts don't occur
  static const int version = 101;

  //first write out the version number
  out.write((char*)&version,sizeof(version));

  static const int nints=5; 
  static const int ndoub=4;
  static const int nlong=0;

  int psf_size = psf.size();
  int spsf_size = spsf.size();
  int ipsf_size = ipsf.size();

  //don't need many variables - most are derived
  int ivars[nints] = {Nkern, Nsub, psf_size, spsf_size, ipsf_size};

  double dvars[ndoub] = {pixscale, sMissingFlux, nMissedPix, psfh};

  long lvars[nlong] = {};

  cout << "variables defined" << endl;

  //write these to file before writing the psfs themselves
  out.write((char*)&nints,sizeof(nints));
  out.write((char*)ivars,sizeof(ivars));

  cout << "ints written" << endl;

  out.write((char*)&ndoub,sizeof(ndoub));
  out.write((char*)dvars,sizeof(dvars));
  cout << "doub written" << endl;

  out.write((char*)&nlong,sizeof(nlong));
  out.write((char*)lvars,sizeof(lvars));
  cout << "long written" << endl;

  //now write out each of the PSFS

  //the standard psf array
  for(int i=0;i<psf_size;i++)
    {
      double tmp = psf[i];
      out.write((char*)&tmp,sizeof(tmp));
    }

  cout << "standard psf written" << endl;

  //the small psf array
  for(int i=0;i<spsf_size;i++)
    {
      double tmp = spsf[i];
      out.write((char*)&tmp,sizeof(tmp));
    }

  cout << "small psf written" << endl;

  //the individual psf array
  for(int i=0;i<ipsf_size;i++)
    {
      double tmp = ipsf[i];
      out.write((char*)&tmp,sizeof(tmp));
    }

  cout << "individual psf written" << endl;

  //We're done - check for errors and then close up

  if(out.fail())
    {
      cerr << __FUNCTION__ << ": Error: There was an unknown error writing to file (" << filename << ")" << endl;
      //don't exit - there may be plenty more that can be done
      status = 2;
    }

  out.close();

  return status;
}

//read in the psfs from a binary file
int PSF::read_psf(string filename)//, string path="")
{
  //cerr << 'testing ' <<filename<<'\n';
  //string psf_path = std::getenv("GULLS_BASE_DIR");

  //string file = psf_path+filename;

  /*should just require filename, need to test*/
  string file = filename;
  int status = 0;
  init = false;

  ifstream in(file.c_str(),ios::in | ios_base::binary);
  if(!in)
    {
      //cerr << 'here2, using ' << filename <<'\n';
      cerr << __FUNCTION__ << ": Error: Could not open psf file for input (" << file << ")" << endl;
      status = 1;
      return status;
    }

  //need to be hard coded by the developer. Any time a change is made increment
  //the version num counter (by hand, here, and in the write function) so that
  // file format conflicts don't occur
  static const int version = 101;

  static const int nints=5; 
  static const int ndoub=4;
  static const int nlong=0;

  int nints_, ndoub_, nlong_;

  int ivars[nints];// = {version, Nkern, sNKern, Nsub, psf_size, spsf_size, 
		   //   ipsf_size};

  double dvars[ndoub];// = {pixscale, sMissingFlux, nMissedPix, psfh};

  long lvars[nlong];// = {};

  int fileversion;
  int psf_size;
  int spsf_size;
  int ipsf_size;

  //start reading in the data - first check the file version
  in.read((char*)&fileversion,sizeof(fileversion));

  if(!in)
    {
      cerr << __FUNCTION__ << ": Error: Problem reading in file (" << filename << ")" << endl;
      status = 2;
      return status;
    }

  if(fileversion!=version)
    {
      cerr << __FUNCTION__ << ": Error: The psf file (" << filename << ") was created using a different version (" << fileversion << ") of the psf class. You will need to recompute the psfs using this version of the psf class - version number (" << version << ")" << endl;
      status = 3;
      return status;
      exit(1);
    }

  //if the version numbers match, assume the file will be valid - don't
  //perform any more error checks until the end

  //read in the descriptive variables
  
  //int ivars[nints];// = {version, Nkern, Nsub, psf_size, spsf_size, 
		   //   ipsf_size};
  in.read((char*)&nints_,sizeof(nints_));
  if(nints_!=nints)
    {
      cerr << __FUNCTION__ << "Wrong number of ints" << endl;
      status = 4;
      return status;
    }
  in.read((char*)ivars,sizeof(ivars));
  int Nkern_ = ivars[0];
  int Nsub_ = ivars[1];
  psf_size = ivars[2];
  spsf_size = ivars[3];
  ipsf_size = ivars[4];

  cout << "ivars: " << Nkern_ << " " << Nsub_ << " " << psf_size << " " << spsf_size << " " << ipsf_size << endl;

  //double dvars[ndoub];// = {pixscale, sMissingFlux, nMissedPix, psfh};
  in.read((char*)&ndoub_,sizeof(ndoub_));
  if(ndoub_!=ndoub)
    {
      cerr << __FUNCTION__ << "Wrong number of doubles" << endl;
      status = 4;
      return status;
    }
  in.read((char*)dvars,sizeof(dvars));
  double pixscale_ = dvars[0];

  //long lvars[nlong];// = {};
  in.read((char*)&nlong_,sizeof(nlong_));
  if(nlong_!=nlong)
    {
      cerr << __FUNCTION__ << "Wrong number of ints" << endl;
      status = 4;
      return status;
    }
  in.read((char*)lvars,sizeof(lvars));

  //perform the required set-up - this function is effectively a constructor
  PSF(Nkern_, Nsub_, pixscale_);

  sMissingFlux = dvars[1];
  nMissedPix = dvars[2];
  psfh = dvars[3];

  cout << "dvars: " << pixscale << " " <<  sMissingFlux << " " <<  nMissedPix << " " <<  psfh << endl;

  //prepare the psf vectors
  psf.resize(psf_size);
  spsf.resize(spsf_size);
  ipsf.resize(ipsf_size);

  //read in the main psf vector
  for(int i=0;i<psf_size;i++)
    {
      double tmp;
      in.read((char*)&tmp,sizeof(tmp));
      psf[i] = tmp;
    }

  //read in the small psf vector
  for(int i=0;i<spsf_size;i++)
    {
      double tmp;
      in.read((char*)&tmp,sizeof(tmp));
      spsf[i] = tmp;
    }

  //read in the individual psf vector
  for(int i=0;i<ipsf_size;i++)
    {
      double tmp;
      in.read((char*)&tmp,sizeof(tmp));
      ipsf[i] = tmp;
    }

  //We're done - check for errors and then close up

  if(in.fail())
    {
      cerr << __FUNCTION__ << ": Error: There was an unknown error reading the file (" << filename << ")" << endl;
      //don't exit - there may be plenty more that can be done
      status = 4;
    }

  in.close();

  init = true; //all present and correct
  return status;
}

int PSF::load_crosstalk(string ifname)
{

  ifstream in(ifname.c_str());
  if(!in)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Could not open crosstalk file " << ifname << endl;
      return 1;
    }

  string line;
  vector<double> data;
  vector<double> ct;

  while(!in.eof())
    {
      getline(in,line);
      split(line,data);

      for(int i=0;i<int(data.size());i++)
	{
	  ct.push_back(data[i]);
	}
    }

  return load_crosstalk(ct);
}

int PSF::load_crosstalk(vector<double> ct)
{
  if(int(ct.size())%int(sqrt(ct.size()))!=0)
    {
      cerr << __FILE__ ": " << __FUNCTION__ << ": Error, crosstalk kernel is not square" << endl;
      return 1;
    }
  
  if(sqrt(int(ct.size())%2==0))
    {
      cerr << __FILE__ ": " << __FUNCTION__ << ": Error, crosstalk kernel has an even number of pixels on a side" << endl;
      return 1;
    }
    
  ctalk = ct;
  return 0;
}

void PSF::apply_crosstalk()
{
  vector<double> tmppsf;

  int soffset;
  int ssoffset;
  int ctoffset = sqrt(ctalk.size())/2; //sqrt(ctalk.size()) will be odd
  int ctside = 2*ctoffset+1;
  int yoffset;
  int syoffset;
  int joffset;
  int jjoffset;
  double psfval;

  int x,y;

  //the main psf
  tmppsf = psf;
  psf.clear();
  psf.resize(psfsize,0);
  spsf.clear();
  spsf.resize(spsfsize,0);

  //kernsize*(Nsub*Sy+Sx) + kernside*(Py+Nkern) + (Px+Nkern);

  //for each subpixel psf
  for(int s=0;s<Nsub*Nsub;s++)
  {
    //soffset = 1;
    soffset = s*kernsize;
    ssoffset = s*skernsize;
      //for each pixel
      for(int j=0;j<kernside;j++)
	{
	  joffset = j*kernside;
	  for(int i=0;i<kernside;i++)
	    {
	      psfval = tmppsf[soffset + joffset + i];

	      //disperse the crosstalk across each of the PSF pixels
	      for(int jj=-ctoffset;jj<=ctoffset;jj++)
		{
		  y = jj+j;
		  if(y>=0 && y<kernside)
		    {
		      yoffset = y*kernside;
		      syoffset = (y-Nkern+sNkern)*skernside;
		      jjoffset = (ctoffset + jj)*ctside;
		      for(int ii=-ctoffset;ii<=ctoffset;ii++)
			{
			  x = ii+i;
			  if(x>=0 && x<kernside)
			    {
			      psf[soffset + yoffset + x] += ctalk[jjoffset+ii+ctoffset] * psfval;

			      if(x-Nkern>=-sNkern && x-Nkern<=sNkern && y-Nkern>=-sNkern && y-Nkern<=sNkern)
				{
				  //add to the smaller kernel as well
				  spsf[ssoffset + syoffset + (x-Nkern+sNkern)] += ctalk[jjoffset+ii+ctoffset] * psfval;
				}
			    } //end if x
			} //end ii
		    } //end if y
		} //end jj
	    } //end i
	} //end j
  } //end s


  //renormalize the PSF
  //for each subpixel psf
  //double sum;
  //for(int s=0;s<Nsub*Nsub;s++)
  //{
  //  sum=0;
  //  soffset = s*kernsize;
  //  for(int i=soffset;i<soffset+kernsize;i++) sum+=psf[i];
  //  for(int i=soffset;i<soffset+kernsize;i++) 
  //    {
  //	psf[i]/=sum;
  //    }
  //}

  /*
  //The small psf
  tmppsf = spsf;
  spsf.clear();
  spsf.resize(tmppsf.size(),0);

  //for each subpixel psf
  for(int s=0;s<Nsub*Nsub;s++)
    {
      soffset = s*skernsize;
      //for each pixel
      for(int j=0;j<skernside;j++)
	{
	  joffset = j*skernside;
	  for(int i=0;i<skernside;i++)
	    {
	      psfval = tmppsf[soffset + joffset + i];

	      //disperse the crosstalk across each of the PSF pixels
	      for(int jj=-ctoffset;jj<=ctoffset;jj++)
		{
		  y = jj+j;
		  if(y>=0 && y<skernside)
		    {
		      yoffset = y*skernside;
		      jjoffset = (ctoffset + jj)*ctside;
		      for(int ii=-ctoffset;ii<=ctoffset;ii++)
			{
			  x = ii+i;
			  if(x>=0 && x<skernside)
			    {
			      spsf[soffset + yoffset + x] += ctalk[jjoffset+ii+ctoffset] * psfval;
			    } //end if x
			} //end ii
		    } //end if y
		} //end jj
	    } //end i
	} //end j
    } //end s

  //renormalize the sPSF
  //for each subpixel spsf
  for(int s=0;s<Nsub*Nsub;s++)
  {
    sum=0;
    soffset = s*skernsize;
    for(int i=soffset;i<soffset+skernsize;i++) sum+=spsf[i];
    for(int i=soffset;i<soffset+skernsize;i++) spsf[i]/=sum;
  }
  */
}
