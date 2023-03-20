#ifndef INTEGERPOWERS
#define INTEGERPOWERS

#include<complex>
#include<iostream>
#include<cstdlib>
#include<cmath>

using namespace std;

template<class T> inline const T sqr(const T a) //one multiply
{
  return a*a;
}

template<class T> inline const T cube(const T a) //two multiply
{
  return a*sqr(a);
}

template<class T> inline const T hCube(const T a) //two multiply
{
  return sqr(sqr(a));
}

template<class T> inline const T iPow(const T a, const int beta)
{
  T b;

  switch(beta)
    {
    case 0:
      return T(1.0);
      break;

    case 1:
      return a;
      break;

    case 2:
      return sqr(a);
      break;

    case 3:
      return cube(a);
      break;

    case 4:
      return hCube(a);
      break;

    case 5:
      return a*sqr(sqr(a));
      break;

    case 6:
      return sqr(cube(a));
      break;

    case 7:
      return sqr(sqr(sqr(a)))/a;
      break;

    case 8:
      return sqr(sqr(sqr(a)));
      break;

    case 9:
      return a*sqr(sqr(sqr(a)));
      break;

    case 10:
      b=sqr(a);
      return b*sqr(sqr(b));
      break;

    case 12:
      return sqr(sqr(cube(a)));
      break;

    case 16:
      return hCube(hCube(a));
      break;

    default:
      if(beta%2==0) return sqr(iPow(a,beta/2));
      else return cube(a)*iPow(a,beta-3);
      //cout << "Error: a^beta, beta=" << beta << " hasn't been implemented" << endl;
      //exit(1);
    }
}

inline const double qAdd(const double a, const double b)
{
  if(a==0&&b==0) return 0;
  return ((abs(a)>=abs(b)) ? abs(a)*sqrt(1.0+sqr(b/a)):abs(b)*sqrt(1.0+sqr(a/b)) );
}

inline const double qAdd(const double a, const double b, const double c)
{
  double aa=abs(a);
  double ab=abs(b);
  double ac=abs(c);

  if(aa>=ab && aa>=ac) return aa*sqrt(1.0+sqr(b/a)+sqr(c/a));
  else if(ab>=aa && ab>=ac) return ab*sqrt(1.0+sqr(a/b)+sqr(c/b));
  else return ac*sqrt(1.0+sqr(a/c)+sqr(b/c));
}

inline void qroot(const double a, const double b, const double c, complex<double>& r1, complex<double>& r2)
{
  double q;
  double disc;
  complex<double> cq;

  if((disc=sqr(b)-4.0*a*c)<0.0)
    {
      //complex roots
      if(b<0.0) cq=-0.5*(b-sqrt(complex<double>(disc,0.0)));
      else cq=-0.5*(b+sqrt(complex<double>(disc,0.0)));

      r1=cq/a;
      r2=c/cq;
    }
  else
    {
      //real roots
      if(b<0.0) q=-0.5*(b-sqrt(disc));
      else q=-0.5*(b+sqrt(disc));

      r1=q/a;
      r2=c/q;
    }

}

inline void qroot(const complex<double> a, const complex<double> b, const complex<double> c, complex<double>& r1, complex<double>& r2)
{
  complex<double> disc;
  complex<double> cq;
  double sign=1.0;

  disc=sqrt(sqr(b)-4.0*a*c);

 if(real(conj(b)*disc)<0.0) sign=-1.0;

 cq = -0.5 * (b + sign*disc);

 r1 = cq/a;
 r2 = c/cq;
}

/*template<class T> inline const T max(const T a_,const T b_)
{
  return (a_>=b_?a_:b_);
}

template<class T> inline const T min(const T a_,const T b_)
{
  return (a_<=b_?a_:b_);
  }*/

#endif /*INTEGERPOWERS*/

