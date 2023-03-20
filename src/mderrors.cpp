#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "split.h"
#include "structures.h"
#include "constants.h"
#include "integerPowers.h"

using namespace std;

vector<double> mderrors(gsl_matrix* cov, int npar, int obs1, int obs2, double alpha, int ffp, struct event* Event, struct slcat *Sources, struct slcat *Lenses)
{
  /*
  ifstream in(filename.c_str());

  if(!in)
    {
      cerr << __FILE__ << ": Error: could not open input file (" << filename << ")" << endl;
      exit(1);
    }

  string line;
  vector<string> data;
  */
  //int stage=0;
  //int npar=0;
  int i=0;
  int j=0;



  gsl_matrix* cij; //cov

  cij = gsl_matrix_calloc(npar,npar); //+1 for color surface brightness scatter
  for(i=0;i<npar;i++)
    {
      for(j=0;j<npar;j++)
	{
	  gsl_matrix_set(cij,i,j,gsl_matrix_get(cov,i,j));
	}
    }

  //set color surface brightness scatter
  //gsl_matrix_set(cij,npar,npar, sqr(0.07));
  

  //Compute the jacobian
  gsl_matrix* jacobian = gsl_matrix_calloc(6,npar);  //initialize to zero

  //Target
  // 0 = Mass
  // 1 = Distance

  //Input
  // 0 = t0      
  // 1 = tE      
  // 2 = u0      
  // 3 = rs      
  // 4 = piEN    
  // 5 = piEE    
  // 6,8,... = F0     
  // 7,9,... = fs

  double M = Lenses->data[Event->lens][MASS] * (ffp ? Event->params[QQ] : 1);
  double piEN = Event->piEN;
  double piEE = Event->piEE;
  double piE = Event->piE;
  double rho = Event->rs;
  double pirel = 1.0/Lenses->data[Event->lens][DIST] - 1.0/Sources->data[Event->source][DIST];

  //cout << M << " " << piEN << " " << piEE << " " << piE << " " << rho << " " << pirel << endl;
  //cout << obs1 << " " << obs2 << endl;

  // dlog(M)/dx

  //The alpha-0.2 are incorrect  -it should be alpha + 0.2 #MP @ 07/11/2021
  //The signs on the pien and piee componants are wrong, and they should be divided by 2 #MP @ 07/11/2021
  //Both these changes not made #MP @ 07/11/2021
  
  gsl_matrix_set(jacobian,0,3, -1.0/rho/ln10          ); //rs
  gsl_matrix_set(jacobian,0,4, piEN/sqr(piE)/ln10 ); //piEN
  gsl_matrix_set(jacobian,0,5, piEE/sqr(piE)/ln10 ); //piEE
  gsl_matrix_set(jacobian,0,6, -2.5*alpha/ln10 ); //F0_1
  gsl_matrix_set(jacobian,0,7, -2.5*alpha/ln10/Event->fs[obs1]); //fs_1
  gsl_matrix_set(jacobian,0,8, 2.5*(alpha-0.2)/ln10  ); //F0_2
  gsl_matrix_set(jacobian,0,9, 2.5*(alpha-0.2)/ln10/Event->fs[obs2]); //fs_2
  //gsl_matrix_set(jacobian,0,npar, ln10       ); //color-surface brightness scatter

  // dlog(pirel)/dx
  gsl_matrix_set(jacobian,1,3, -1.0/rho/ln10          ); //rs
  gsl_matrix_set(jacobian,1,4, -piEN/sqr(piE)/ln10 ); //piEN
  gsl_matrix_set(jacobian,1,5, -piEE/sqr(piE)/ln10 ); //piEE
  gsl_matrix_set(jacobian,1,6, -2.5*alpha/ln10 ); //F0_1
  gsl_matrix_set(jacobian,1,7, -2.5*alpha/ln10/Event->fs[obs1]); //fs_1
  gsl_matrix_set(jacobian,1,8, 2.5*(alpha-0.2)/ln10  ); //F0_2
  gsl_matrix_set(jacobian,1,9, 2.5*(alpha-0.2)/ln10/Event->fs[obs2]); //fs_2
  //gsl_matrix_set(jacobian,1,npar, ln10      ); //color-surface brightness scatter

  //dV-I/dx
  gsl_matrix_set(jacobian,2,6, 2.5/ln10); //F0_1
  gsl_matrix_set(jacobian,2,7, 2.5/ln10/Event->fs[obs1]); //fs_1
  gsl_matrix_set(jacobian,2,8, -2.5/ln10); //F0_2
  gsl_matrix_set(jacobian,2,9, -2.5/ln10/Event->fs[obs2]); //fs_2

  //dthetaE/dx
  gsl_matrix_set(jacobian,3,3, -1.0/rho/ln10); //rs
  gsl_matrix_set(jacobian,3,6, -2.5*alpha/ln10 ); //F0_1
  gsl_matrix_set(jacobian,3,7, -2.5*alpha/ln10/Event->fs[obs1]); //fs_1
  gsl_matrix_set(jacobian,3,8, 2.5*(alpha-0.2)/ln10  ); //F0_2
  gsl_matrix_set(jacobian,3,9, 2.5*(alpha-0.2)/ln10/Event->fs[obs2]); //fs_2

  //dpiE/dx
  gsl_matrix_set(jacobian,4,4, piEN/sqr(piE)/ln10 ); //piEN
  gsl_matrix_set(jacobian,4,5, piEE/sqr(piE)/ln10 ); //piEE


  //dphipiE/dx
  gsl_matrix_set(jacobian,5,4, -piEE/sqr(piE) * 180.0/pi ); //piEN
  gsl_matrix_set(jacobian,5,5, piEN/sqr(piE) * 180.0/pi ); //piEE




  //compute the covariance matrix

  gsl_matrix* temp = gsl_matrix_calloc(npar,6);
  gsl_matrix* Mpirel = gsl_matrix_calloc(6,6);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, cij, jacobian, 0.0, temp);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, jacobian, temp, 0.0, Mpirel);



  vector<double> retval(7);
  //M (M-pirel cov) pirel, color, thetaE, |piE|, arg(piE)
  retval[0] = gsl_matrix_get(Mpirel,0,0);
  retval[1] = gsl_matrix_get(Mpirel,0,1);
  retval[2] = gsl_matrix_get(Mpirel,1,1);
  retval[3] = gsl_matrix_get(Mpirel,2,2);
  retval[4] = gsl_matrix_get(Mpirel,3,3);
  retval[5] = gsl_matrix_get(Mpirel,4,4);
  retval[6] = gsl_matrix_get(Mpirel,5,5);


  //clean-up
  gsl_matrix_free(cij);
  gsl_matrix_free(jacobian);
  gsl_matrix_free(temp);
  gsl_matrix_free(Mpirel);

  return retval;
}


vector<double> piEdirerrors(gsl_matrix* cov, int npar, struct event* Event, struct slcat *Sources, struct slcat *Lenses)
{
  /*
  ifstream in(filename.c_str());

  if(!in)
    {
      cerr << __FILE__ << ": Error: could not open input file (" << filename << ")" << endl;
      exit(1);
    }

  string line;
  vector<string> data;
  */
  //int stage=0;
  //int npar=0;
  int i=0;
  int j=0;



  gsl_matrix* cij; //cov

  cij = gsl_matrix_calloc(npar,npar); //+1 for color surface brightness scatter
  for(i=0;i<npar;i++)
    {
      for(j=0;j<npar;j++)
	{
	  gsl_matrix_set(cij,i,j,gsl_matrix_get(cov,i,j));
	}
    }

  //Target
  // 0 = log(piE)
  // 1 = theta_piE

  //Input
  // 0 = t0      
  // 1 = tE      
  // 2 = u0      
  // 3 = rs      
  // 4 = piEN    
  // 5 = piEE    
  // 6,8,... = F0     
  // 7,9,... = fs

  double piEN = Event->piEN;
  double piEE = Event->piEE;
  double piE = Event->piE;
  double rho = Event->rs;

  //Compute the jacobian
  gsl_matrix* jacobian = gsl_matrix_calloc(2,npar);  //initialize to zero

  gsl_matrix_set(jacobian,0,4,piEN/sqr(piE)/ln10);
  gsl_matrix_set(jacobian,0,5,piEE/sqr(piE)/ln10);
  gsl_matrix_set(jacobian,1,4,-piEE/sqr(piE));
  gsl_matrix_set(jacobian,1,5,piEN/sqr(piE));


  //compute the covariance matrix
  gsl_matrix* temp = gsl_matrix_calloc(npar,2);
  gsl_matrix* piEdir = gsl_matrix_calloc(2,2);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, cij, jacobian, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, jacobian, temp, 0.0, piEdir);

  vector<double> retval(3);
  retval[0] = gsl_matrix_get(piEdir,0,0);
  retval[1] = gsl_matrix_get(piEdir,0,1)*180.0/pi;
  retval[2] = gsl_matrix_get(piEdir,1,1)*sqr(180.0/pi);

  //clean-up
  gsl_matrix_free(cij);
  gsl_matrix_free(jacobian);
  gsl_matrix_free(temp);
  gsl_matrix_free(piEdir);

  return retval;
}
