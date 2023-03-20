#ifndef SINGLELENS
#define SINGLELENS

#include "integerPowers.h"

inline double pacAmp(double r)
{
  return (r*r+2.0)/(r*qAdd(r,2.0));
}

inline double pacAmp(double t, double t0, double tE, double u0)
{
  return pacAmp(qAdd(u0,(t-t0)/tE));
}

inline double pacAmp(complex<double> rs)
{
  return pacAmp(abs(rs));
}

inline double pacMag(double t, double t0, double tE, double u0, double m0, double fs=1.0)
{
  return m0 - 2.5*log10(fs*pacAmp(t,t0,tE,u0)+(1.0-fs));
}

inline double pacMag(complex<double> rs, double m0, double fs=1.0)
{
  return m0 - 2.5*log10(fs*pacAmp(rs)+(1.0-fs));
}

inline double pacMag(double r, double m0, double fs=1.0)
{
  return m0-2.5*log10(fs*pacAmp(r)+(1.0-fs));
}

inline double pspl_dmudu(double u)
{
  double v = u*u;
  double a = v+2.0;
  double b = sqrt(v*(v+4));
  return 2.0*u*(1.0-sqr(a/b))/b;
}

#endif /* SINGLELENS */
