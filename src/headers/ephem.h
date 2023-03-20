//Class for calculating orbits of planets etc.
//Implementation of E. M. Standish & J. G. Williams, Orbital Ephemerides 
//of the Sun, Moon and Planets, in: Explanatory Supplement to the 
//Astronomical Almanac, Chapter 8; http://iau-comm4.jpl.nasa.gov/XSChap8.pdf
//ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/ExplSupplChap8.pdf

#ifndef EPHEMERIS_HEADER

#include<vector>
#include<cmath>

#include "constants.h"

using namespace std;

class orbitalElements 
{
 private:

  //numerical accuracy required for solution of Kepler's equation
  static const double tol = 1.0e-6*3.1415926535897932384626433/180.0; 

  //orbital elements and their time derivatives (in per year)
  //present value, value at epoch, rate of change
  //values will be stored in radians, regardless of input
  double epoch;
  double a, a0, da; //semimajor axis           (AU)
  double e, e0, de; //eccentricity             (rad)
  double I, I0, dI; //inclination              (deg)
  double L, L0, dL; //mean longitude           (deg)
  double w, w0, dw; //longitude of perihelion  (deg) \varomega
  double O, O0, dO; //longitude of asc node    (deg) \Omega
  double b, c, s, f; //aditional terms for outer planet long range ephem (deg)

  double t; //current time in days
  double T; //current time in centuries

  double W; //argument of perihelion \omega
  double M; //mean anomaly
  double E; //eccentric anomaly

  void storerad()
  {
    I0 *= d2r; dI *= d2r;
    L0 *= d2r; dL *= d2r;
    w0 *= d2r; dw *= d2r;
    O0 *= d2r; dO *= d2r;
    b *= d2r;  c  *= d2r; s *= d2r; f *= d2r;
  };

  //functions for computing positions
  void meananomaly();
  void eccentricanomaly();
  void heliocoords();
  void eclcoords();
  void eqcoords();

  //utility functions

  //  double dot(vector<double> x, vector<double> y)
  //{
  // return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  //}

 public:

  static const double J2000 = 2451545.0;
  static const double cosob = 0.9174821392;
  static const double sinob = 0.397776978;
  static const double AUday2kms=1731.45684;
  static const double day2s = 86400.0;
  static const double aum = 1.49597871e11; //AU in m
  static const double aukm = 1.49597871e8; //AU in km
  static const double sindngp=0.455983776;
  static const double cosdngp=0.889988087;
  static const double sinangp=-0.222560701;
  static const double cosangp=0.974918834;
  static const double rangp=3.36603292;
  static const double ra0=4.936838322;
  static const double ll0=0.574736922;
 
  double xh, yh, zh; //heliocentric coordinates in the plane of orbit
  double xecl, yecl, zecl; //coordinates in the J2000 ecliptic plane
  double xeq, yeq, zeq; //coordinates in J2000 equatorial plane (ICRF)

  //constructors
  orbitalElements()
    {
      //default to earth-moon barycenter
      earthmoonbary();
    };

  //if you don't care about the long term
  orbitalElements(double a, double e, double I, double L, double w, double O, double dL_, double epoch_=J2000)
    {
      epoch=epoch_;
      a0=a; da=0;
      e0=e; de=0;
      I0=I; dI=0;
      L0=L; dL=dL_;
      w0=w; dw=0;
      O0=O; dO=0;
      b=0; c=0; s=0; f=0;

      storerad();
    };

  //long term, outer planets
  orbitalElements(vector<double> elements, vector<double> derivatives, double epoch_=J2000, double bb=0, double cc=0, double ss=0, double ff=0)
    {
      epoch=epoch_;
      a0=elements[0]; da=derivatives[0];
      e0=elements[1]; de=derivatives[1];
      I0=elements[2]; dI=derivatives[2];
      L0=elements[3]; dL=derivatives[3];
      w0=elements[4]; dw=derivatives[4];
      O0=elements[5]; dO=derivatives[5];
      b=bb; c=cc; s=ss; f=ff;

      storerad();
    };

  //long term, outer planets
  orbitalElements(double a00, double e00, double I00, double L00, double w00, 
		  double O00,
		  double ad, double ed, double Id, double Ld, double wd, 
		  double Od,
		  double epoch_=J2000,
		  double bb=0, double cc=0, double ss=0, double ff=0)
    {
      epoch = epoch_;
      a0=a00; da=ad;
      e0=e00; de=ed;
      I0=I00; dI=Id;
      L0=L00; dL=Ld;
      w0=w00; dw=wd;
      O0=O00; dO=Od;
      b=bb; c=cc; s=ss; f=ff;

      storerad();
    };

  //common planets/orbits

  //Earth-moon barycenter
  void earthmoonbary();

  //Earth about the Earth-moon barycenter - epoch: Jan 1 2022
  void earth();

  //L2 orbit
  void earthl2();

  //Stationary orbit
  void blank();

  //suitable for long term calculations
  void mars();

  //suitable for short term calculations
  void jupiter();

  //jwst orbit
  void jwst(double phase=0);

  //geosynchronous and geostationary Earth orbits
  void geosynch(double inc=0, double phase=0); //inc is referenced to equator
  void geostat(double phase=0)
  {
    geosynch(0,phase);
  }

  //Approximations of a lissajous orbit
  void lissajous(double phase=0);
  void lissajousxy(double phase=0);
  void lissajousz(double phase=0);

  //destuctor
  ~orbitalElements(){};

  //functions

  void settime(double jd);

  void compute();

  //for returning computed positions

  //heliocentric in the orbital plane of the planet
  void heliocentric(vector<double>* helio)
  {
    helio->resize(3);
    (*helio)[0] = xh;
    (*helio)[1] = yh;
    (*helio)[2] = zh;    
  };

  //heliocentric in the ecliptic plane
  void ecliptic(vector<double>* ecl)
  {
    ecl->resize(3);
    (*ecl)[0] = xecl;
    (*ecl)[1] = yecl;
    (*ecl)[2] = zecl;    
  };

  //heliocentric in the equatorial plane
  void j2000(vector<double>* jxyz)
  {
    jxyz->resize(3);
    (*jxyz)[0] = xeq;
    (*jxyz)[1] = yeq;
    (*jxyz)[2] = zeq;    
  };

  //Compute this orbit and perturb it with other orbits
  void perturb(double t, vector<orbitalElements*> perturbers);
  
  //perturb this orbit with other orbits
  void perturb(orbitalElements* perturbation);

  //view the orbit from a star at ra,dec
  //coordinates of output are N,E,Distance(AU)
  void viewfrom(double ra, double dec, vector<double>* view);
  void viewfrom_old(double ra, double dec, vector<double>* view);

  //compute and view the orbit from a star at ra,dec
  void viewfrom(double t, double ra, double dec, vector<double>* view);

  //velocity
  //coordinates of output are vN,vE(AU/d),vrad(km/s)
  void velocity(double t, double ra, double dec, vector<double>* v, double tstep=1e-6);

  //acceleration
  //coordinates of output are aN,aE(AU/d^2),arad(km/s^2)
  void acceleration(double t, double ra, double dec, vector<double>* acc, double tstep=1e-6);

  //change in proper motion - subtract this from the heliocentric pm
  void eqpropermotion(double t, double ra, double dec, double dist, vector<double>* pm, double tstep=0.001);

  //change in proper motion - subtract this from the heliocentric pm
  void galpropermotion(double t, double ll, double bb, double dist, vector<double>* pm, double tstep=0.001);

  void eq2gal(vector<double> eq, vector<double>* gal);
  void gal2eq(vector<double> gal, vector<double>* eq);
  void pmeq2gal(vector<double> eq, vector<double> pmeq, vector<double>* pmgal);
  void pmgal2eq(vector<double> gal, vector<double> pmgal, vector<double>* pmeq);

  double get_epoch()
  {
    return epoch;
  };
  
};

#define EPHEMERIS_HEADER
#endif
