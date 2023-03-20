#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

#include "rbf.h"
#include "integerPowers.h"
#include "split.h"

using namespace std;

void rbf::construct_interpolation(vector<double>* _x, vector<double>* _y, double _c)
{
 
  if(constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to construct an interpolation that is already constructed. Use clear before reconstructing." << endl;
      exit(1);
    }
     
  ndata = _y->size();
  if(_x->size() % ndata > 0)
    {
      cerr << __FILE__ << ": Error: Dimensions of abscissa and value vectors do not match: " << _x->size() << "," << _y->size() << endl;
      exit(1);
    }
      
  c=_c;
  c2=c*c;
  ndim = _x->size()/ndata;

  //cout << "ndim = " << ndim << endl;

  gsl_vector* y = gsl_vector_alloc(ndata);
  lambda = gsl_vector_alloc(ndata);
  constructed=1;
  scale.resize(ndim);
  x.resize(ndim*ndata);

  //load the data
  for(int i=0;i<ndata;i++)
    {
      gsl_vector_set(y,i,(*_y)[i]);
      gsl_vector_set(lambda,i,(*_y)[i]);

      int ii=i*ndim;

      for(int k=0;k<ndim;k++)
	{
	  x[ii+k] = (*_x)[ii+k];
	}
    }

  //Scale each axis
      
  for(int k=0;k<ndim;k++) scale[k]=0;
  for(int i=0;i<ndata;i++)
    {
      int ii=ndim*i;
      for(int k=0;k<ndim;k++) scale[k] += x[ii+k];
    }
  for(int k=0;k<ndim;k++) scale[k]/=double(ndata);

  for(int i=0;i<ndata;i++)
    {
      int ii=ndim*i;
      for(int k=0;k<ndim;k++) 
	{
	  x[ii+k]/=scale[k];
	  //cout << x[ii+k] << endl;
	}
    }

  //Now work out a sensible value for c
  if(c==RBF_DEFAULT_C)
    {
      double cguess=0;
      
      for(int i=0;i<ndata;i++)
	{
	  int ii=i*ndim;
	  double mindist=1e50;
	  for(int j=0;j<ndata;j++)
	    {
	      double sum=0;
	      int jj=j*ndim;
	      if(j!=i)
		{
		  for(int k=0;k<ndim;k++)
		    {
		      sum+=sqr(x[ii+k]-x[jj+k]);
		    }
		  if(sum<sqr(mindist)) mindist=sqrt(sum);
		}
	    }
      
	  cguess+=mindist;
	}

      c = cguess/double(ndata);
      //cout << "#The best guess value for c is " << c << endl;
    }
	         

  //calculate the matrix of the norms
  gsl_matrix* Aij = gsl_matrix_alloc(ndata,ndata);
  gsl_matrix* LU = gsl_matrix_alloc(ndata,ndata);

  for(int i=0;i<ndata;i++)
    {
      int ii=ndim*i;
      for(int j=0;j<ndata;j++)
	{
	  if(i>j) gsl_matrix_set(Aij,i,j,gsl_matrix_get(Aij,j,i));
	  else
	    {
	      int jj=ndim*j;
	      double sum=1.0;
	      for(int k=0;k<ndim;k++)
		{
		  sum+=sqr((x[ii+k]-x[jj+k])/c);
		}
	      gsl_matrix_set(Aij,i,j,c*sqrt(sum));
	    }
	  gsl_matrix_set(LU,i,j,gsl_matrix_get(Aij,i,j));
	}
    }

  for(int i=0;i<ndata;i++)
    {
      for(int j=0;j<ndata;j++)
	{
	  //cout << gsl_matrix_get(LU,i,j) << " ";
	}
      //cout << endl;
    }
  //cout << endl << endl;

  //calculate the LU decomposition
  gsl_permutation* permutation = gsl_permutation_alloc(ndata);
  int signum;
      
  gsl_linalg_LU_decomp(LU, permutation, &signum);

  for(int i=0;i<ndata;i++)
    {
      for(int j=0;j<ndata;j++)
	{
	  //cout << gsl_matrix_get(LU,i,j) << " ";
	}
      //cout << endl;
    }
  //cout << endl << endl;

      
  //calculate the coefficients of the interpolation
  gsl_linalg_LU_solve(LU, permutation, y, lambda);

  //refine the solution
  //gsl_vector* residual = gsl_vector_alloc(ndata);
  
  // for(int r=0;r<3;r++)
  //   {
  //     gsl_linalg_LU_refine(Aij,LU,permutation,y,lambda,residual);
  //     cout << "residual(" << r << "): ";
  //     for(int i=0;i<ndata;i++) cout << gsl_vector_get(residual,i) << " ";
  //     cout << endl;
  //   }

  // //check for a large residual
  // int largeResidual=0;
  // for(int i=0;i<ndata;i++)
  //   {
  //     if(abs(gsl_vector_get(residual,i))>1.0e-5) 
  // 	cerr << "Warning: Large residual in interpolation element " << i << ": " << gsl_vector_get(residual,i) << endl;
  //   }

  //delete anything we don't need
  gsl_matrix_free(Aij);
  gsl_matrix_free(LU);
  gsl_vector_free(y);
  gsl_permutation_free(permutation);
  //gsl_vector_free(residual);

  //cout << "ndim = " << ndim << endl;
}

double rbf::interpolate(double _x[], double _y[], int _ndata, int _ndim, double _c)
{
  if(constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to construct an interpolation that is already constructed. Use clear before reconstructing." << endl;
      exit(1);
    }

  vector<double> __x;
  vector<double> __y;

  for(int i=0;i<_ndata;i++)
    {
      int ii=_ndim*i;
      __y.resize(_ndata);
      __x.resize(_ndim*_ndata);
	  
      __y[i]=_y[i];
      for(int k=0;k<_ndim;k++)
	{
	  __x[ii+k] = _x[ii+k];
	}
    }

  //cout << "Alternative constructor" << endl;
  constructed=0;
  construct_interpolation(&__x, &__y, _c);

  return c;
  //cout << "ndim = " << ndim << endl;
      
}

double rbf::evaluate(vector<double>* _xpoint)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  double sum=0;
  double rsum;

  if(int(_xpoint->size())!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _xpoint->size() = " << _xpoint->size() << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  vector<double> xpoint;

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back((*_xpoint)[k]/scale[k]);
    }

  for(int i=0;i<ndata;i++)
    {
      rsum=1;
      int ii=i*ndim;

      for(int k=0;k<ndim;k++)
	{
	  rsum+=sqr((xpoint[k]-x[ii+k])/c);
	}

      sum += gsl_vector_get(lambda,i)*c*sqrt(rsum);
    }

  return sum;
}
 
double rbf::evaluate(double _xpoint[], int _ndim)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  vector<double> xpoint;

  if(_ndim!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _ndim = " << _ndim << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back(_xpoint[k]);
    }

  return evaluate(&xpoint);
}

double rbf::derivative(vector<double>* _xpoint, int dim)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  double sum=0;
  double rsum;

  if(int(_xpoint->size())!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _xpoint->size() = " << _xpoint->size() << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  vector<double> xpoint;

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back((*_xpoint)[k]/scale[k]);
    }

  for(int i=0;i<ndata;i++)
    {
      rsum=1;
      int ii=i*ndim;

      for(int k=0;k<ndim;k++)
	{
	  rsum+=sqr((xpoint[k]-x[ii+k])/c);
	}

      sum += gsl_vector_get(lambda,i)*(xpoint[dim]-x[ii+dim])/(c*sqrt(rsum));
    }

  return sum/scale[dim];
}
 
double rbf::derivative(double _xpoint[], int dim, int _ndim)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  vector<double> xpoint;

  if(_ndim!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _ndim = " << _ndim << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back(_xpoint[k]);
    }

  return derivative(&xpoint, dim);
}

complex<double> rbf::derivative_fdf(vector<double>* _xpoint, int dim)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  double fsum=0;
  double dfsum=0;
  double rsum;

  if(int(_xpoint->size())!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _xpoint->size() = " << _xpoint->size() << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  vector<double> xpoint;

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back((*_xpoint)[k]/scale[k]);
    }

  for(int i=0;i<ndata;i++)
    {
      rsum=1;
      int ii=i*ndim;
      double l = gsl_vector_get(lambda,i);

      for(int k=0;k<ndim;k++)
	{
	  rsum+=sqr((xpoint[k]-x[ii+k])/c);
	}

      fsum += l*c*sqrt(rsum);
      dfsum += l*(xpoint[dim]-x[ii+dim])/(c*sqrt(rsum));
    }

  return complex<double>(fsum,dfsum/scale[dim]);
}
 
complex<double> rbf::derivative_fdf(double _xpoint[], int dim, int _ndim)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to evaluate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  vector<double> xpoint;

  if(_ndim!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _ndim = " << _ndim << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  for(int k=0;k<ndim;k++)
    {
      xpoint.push_back(_xpoint[k]);
    }

  return derivative_fdf(&xpoint, dim);
}

double rbf::integrate(vector<double>* xmin, vector<double>* xmax)
{
  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to integrate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  if(int(xmin->size())!=ndim || int(xmax->size())!=ndim)
    {
      cerr << __FILE__ << ": Error: Integration limit dimensions do not match the input dimensions" << endl;
    }

  //rescale the inputs

  for(int k=0;k<ndim;k++)
    {
      (*xmin)[k] /= scale[k];
      (*xmax)[k] /= scale[k];
    }

  if(ndim == 1)
    {
      //result is calculated...
      double sum=0;
      double dx;
      double root;

      for(int i=0;i<ndata;i++)
	{
	  //add the xmax part
	  dx = (*xmax)[0] - x[i];
	  root = c*sqrt(1.0+sqr(dx/c));
	  sum += gsl_vector_get(lambda,i) * (dx*root + c2*log(root+dx));

	  //subtract the xmin part
	  dx = (*xmin)[0] - x[i];
	  root = c*sqrt(1.0+sqr(dx/c));
	  sum -= gsl_vector_get(lambda,i) * (dx*root + c2*log(root+dx));
	}

      sum *= 0.5*scale[0];

      return sum;
    }
  else
    {
      //can only do 1d integrals - use numerics after that

      //look for the one dimension which has size > 0

      vector<double> size;
      int ndimensions=0;
      int intdim=0;
	
      for(int k=0;k<ndim;k++)
	{
	  if(abs((*xmax)[k]-(*xmin)[k])>1e-10)
	    {
	      ndimensions++;
	      intdim=k;
	    }
	}

      if(ndimensions!=1)
	{
	  cerr << __FILE__ << ": Error: There are " << ndimensions << " integration dimensions. Can only have one. Use a numerical integrator for dimensions greater than one." << endl;
	  exit(1);
	}

      //now calculate the integral

      double sum=0;
      double dx;
      double root;
      double ceff, ceff2;

      for(int i=0;i<ndata;i++)
	{
	  int ii=ndim*i;

	  //calculate the effective c
	  ceff2=c2;
	  for(int k=0;k<ndim;k++)
	    {
	      if(k!=intdim) ceff2+=sqr((*xmin)[k]-x[ii+k]);
	    }
	  ceff = sqrt(ceff2);

	  //add the xmax part
	  dx = (*xmax)[intdim] - x[ii+intdim];
	  root = ceff*sqrt(1.0+sqr(dx/ceff));
	  sum +=  gsl_vector_get(lambda,i) * (dx*root + ceff2*log(root+dx));

	  //subtract the xmin part
	  dx = (*xmin)[intdim] - x[ii+intdim];
	  root = ceff*sqrt(1.0+sqr(dx/ceff));
	  sum +=  gsl_vector_get(lambda,i) * (dx*root + ceff2*log(root+dx));
	}

      sum *= 0.5*scale[intdim];

      return sum;
	
    }
}

double rbf::integrate(double xmin[], double xmax[], int _ndim)
{

  if(!constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to integrate from an interpolation that has not been constructed" << endl;
      exit(1);
    }

  if(_ndim!=ndim)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: _ndim = " << _ndim << " does not match the ndim that was originally provided: ndim = " << ndim << endl;
      exit(1);
    }

  vector<double> _xmin;
  vector<double> _xmax;

  for(int k=0;k<_ndim;k++)
    {
      _xmin.push_back(xmin[k]);
      _xmax.push_back(xmax[k]);
    }

  return integrate(&_xmin, &_xmax);
}

//Load data table from file
double rbf::interpolate(string input, double _c)
{

  if(constructed)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to construct an interpolation that is already constructed. Use clear before reconstructing." << endl;
      exit(1);
    }

  ifstream in(input.c_str());
  if(!in)
    {
      cerr << "Could not open input file: " << input << endl;
      exit(1);
    }

  vector<double> x_;
  vector<double> y_;

  vector<string> data;
  string line;

  int nfields=-1;


  //load in the data
  while(!in.eof())
    {
      getline(in,line);
      if(line.find_first_of("#ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz")==string::npos)
	{
	  split(line,data);
	  if(nfields==-1&&int(data.size())>1) nfields = data.size();
	  if(int(data.size())==nfields)
	    {
	      for(int i=0;i<int(data.size()-1);i++)
		{
		  x_.push_back(atof(data[i].c_str()));
		}
	      y_.push_back(atof(data[data.size()-1].c_str()));
	    }
	}
    }

  //compute the interpolation

  return interpolate(&x_,&y_,_c);
}
