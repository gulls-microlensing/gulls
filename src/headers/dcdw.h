#ifndef DCDW_HEADER
#define DCDW_HEADER

#include<cmath>

#include "integerPowers.h"

inline double drw(double q)
{
  return sqrt(pow((1.0+pow(q,1.0/3.0)),3.0)/(1.0+q));
}
  
inline double dcr(double q)
{
  double y, a_, g;
  
  a_ = 27.0*q/sqr(1.0+q);
  
  g = double(pow( double(2.0/(a_*(-27.0+sqrt(729.0-108.0*a_)+a_*(18.0-2.0*a_)))) ,(1.0/3.0)));
  
  y = (1.0/3.0)*(3.0-a_-g*a_*(6.0-a_)+1.0/g);
  
  return pow(y,0.25);
  
}

inline int causticTopology(double s,double q)
{
  //close = 1
  //resonant = 2
  //wide = 3
  //  return (d_<dcr(q)? 1 : (d_<drw(q) ? 2 : 3) );
  return (s>drw(q) ? 3 : (s>dcr(q) ? 2 : 1) );
}

#endif /*DCDW_HEADER*/
