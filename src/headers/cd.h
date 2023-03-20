#ifndef CDTYPECAST
#define CDTYPECAST

#include<complex>

#include "integerPowers.h"

using namespace std;

typedef complex<double> cd;

inline const cd zzb(const cd& z_)
{
  return sqr(real(z_)) - sqr(imag(z_));
}

struct cdless {
  bool operator() (const cd& lhs, const cd& rhs) const
  {
    if(real(lhs)<real(rhs)) return true;
    else if(real(lhs)==real(rhs) && imag(lhs)<imag(rhs)) return true;
    else return false;
  }
};

#endif /*CDTYPECAST*/
