#ifndef BINARYMAG_HEADER
#define BINARYMAG_HEADER

///////////////////////////////////////////////////////////////
//
//             WARNING!!!!!
//
//               Does not use centre of mass frame
//               Uses frame where m1 z2 + m2 z1 = 0
//               Use comframe to convert source position
//
///////////////////////////////////////////////////////////////

#include "cd.h"

struct lensParameters
{
  double m1, m2, z1, z2;

  lensParameters(double m1_, double m2_, double z1_, double z2_)
  {
    m1=m1_; m2=m2_; z1=z1_; z2=z2_;

    if(abs(m1+m2-1.0)>1.0e-10||abs(m1*z2+m2*z1)>1.0e-10)
      {
	cerr << "lensParameters: Error: Total mass does not equal 1.0" << endl;
	exit(1);
      }
  };

  lensParameters(double d=0.5, double q=0.5)
  {
    m1 = 1.0/(1.0+q);
    m2 = q*m1;

    z1 = -m1*d;
    z2 = m2*d;
  };
};

inline cd comframe(cd zscom,lensParameters* l)
{
  //converts from centre of mass frame to working frame
  return zscom - (l->m1-l->m2)*(l->z2-l->z1);
}

inline cd emarfmoc(cd zswork,lensParameters* l)
{
  //converts from working frame to centre of mass frame
  return zswork + (l->m1-l->m2)*(l->z2-l->z1);
}

void coefficients(cd zs, lensParameters* l, cd coeffs[6]);

double magnification(cd zs, int& nimgs, cd imgs[5], int& flag, lensParameters* l);

void lensEquation(cd zs, cd& z, cd& df, double& jac, cd& dz,lensParameters* l);

bool analyticCaustic(double phi, lensParameters* l, cd cc[4], cd ca[4]);
bool newtonCaustic(double phi, lensParameters* l, cd cc[4], cd ca[4]);

inline bool checkCriticalCurve(cd z, lensParameters* l)
{
  double check=abs(l->m1/sqr(z-l->z1) + l->m2/sqr(z-l->z2));
  if(abs(check-1.0)>1.0e-7)
    {
      return true;
    }
  else return false;
}

inline cd lensEquation(cd z, lensParameters* l)
{
  return z - l->m1/conj(z-l->z1) - l->m2/conj(z-l->z2);
}

bool lageurCaustic(double phi, lensParameters * l, cd cc[4], cd ca[4]);

int lageurCusp(double lambda, lensParameters * l, cd cusp[6]);

#endif /*BINARYMAG_HEADER*/
