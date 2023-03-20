#ifndef MABULS_CROIN_HEADER

using namespace std;

#include "constants.h"
#include "dcdw.h"

void usecroin(double s, double q, double tc, double tE, double uc, double alpha, double* t0, double* u0, double* ucmax=NULL, double* xc=NULL, double* yc=NULL)
{
  double rcroin;
  double xcroin;
  double ycroin;

  double arad = alpha*pi/180;
  double ca = cos(arad);
  double sa = sin(arad);

  if(causticTopology(s,q)==2)
    {
      xcroin=0;
      ycroin=0;
      rcroin=4.5*pow(q,0.25);
      //rcroin=pow(q,0.25);
    }
  else
    {
      if(s<0.1)
	{
	  xcroin = (s - 1.0/s)*(1-q)/(1+q);
	  ycroin = sqrt(q)/(1+q)*(2.0/s - s);
	  rcroin=1e-10;
	}
      else if(s<1)
	{
	  xcroin = (s - 1.0/s)*(1-q)/(1+q);
	  ycroin=0;
	  rcroin = 2*sqrt(q)/(s*sqrt(1+s*s));
	}
      else
	{
	  xcroin = (s - 1.0/(s*(1+q))) - q*s/(1+q);
	  ycroin = 0;
	  rcroin = 2*sqrt(q);
	}
      xcroin += s*q/(1+q);
      rcroin *= 2 + min(45*s*s,80/s/s);
      //xcroin += s/sqr(1+q) - s*q/(1+q);
      
    }

  *t0 = tc - tE*(xcroin*ca + ycroin*sa);
  *u0 = uc*rcroin - (xcroin*sa - ycroin*ca);
  if(ucmax!=NULL) *ucmax = rcroin;
  if(xc!=NULL) *xc = xcroin;
  if(yc!=NULL) *yc = ycroin;
}

void croinparam(double s, double q, double t0, double tE, double u0, double alpha, double* tc, double* uc, double* ucmax=NULL, double* xc=NULL, double* yc=NULL)
{
  double rcroin;
  double xcroin;
  double ycroin;

  double arad = alpha*pi/180;
  double ca = cos(arad);
  double sa = sin(arad);

  if(causticTopology(s,q)==2)
    {
      xcroin=0;
      ycroin=0;
      rcroin=4.5*pow(q,0.25);
      //rcroin=pow(q,0.25);
    }
  else
    {
      if(s<0.1)
	{
	  xcroin = (s - 1.0/s)*(1-q)/(1+q);
	  ycroin = sqrt(q)/(1+q)*(2.0/s - s);
	  rcroin=1e-10;
	}
      else if(s<1)
	{
	  xcroin = (s - 1.0/s)*(1-q)/(1+q);
	  ycroin=0;
	  rcroin = 2*sqrt(q)/(s*sqrt(1+s*s));
	}
      else
	{
	  xcroin = (s - 1.0/(s*(1+q))) - q*s/(1+q);
	  ycroin = 0;
	  rcroin = 2*sqrt(q);
	}
      xcroin += s*q/(1+q);
      rcroin *= 2 + min(45*s*s,80/s/s);
      //xcroin += s/sqr(1+q) - s*q/(1+q);
      
    }

  *tc = t0 + tE*(xcroin*ca + ycroin*sa);
  *uc = (u0 + (xcroin*sa - ycroin*ca))/rcroin;
  if(ucmax!=NULL) *ucmax = rcroin;
  if(xc!=NULL) *xc = xcroin;
  if(yc!=NULL) *yc = ycroin;
}

#define MABULS_CROIN_HEADER
#endif
