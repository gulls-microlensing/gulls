#ifndef RANDOM_HEADER

#include<math.h>

float ran1(long *idum);
float gasdev(long *idum);
float ran3(long *idum);
float ran2(long *idum);
double gammln(double xx);
double poisson(double mean, long* idum);

inline double normal(double mean, double sigma, long* idum)
{
  return mean + sigma*gasdev(idum);
}

inline double logNormal(double mean, double sigma, long* idum)
{
  return exp(log(mean)+sigma*gasdev(idum));
}

inline long randint(long min, long max, long* idum)
{
  //inclusively between
  return min + int(ran2(idum)*(max-min+1));
}

#define RANDOM_HEADER
#endif
