#ifndef ZODIACAL_LIGHT_HEADER

#include<vector>
#include<iostream>
#include<cmath>

#include "rbf.h"

using namespace std;

struct felong_const{
  //f(elong) = a02 + b02 g(elong) + a25 + b25 g(elong)
  //where a = int_lmin^lmax f_filter(l) f_abs(l) dl
  //and b = int_lmin^lmax f_filter(l) f_abs(l) log(l) dl
  //and g(elong) has a simple form
  
  double a02;
  double b02;
  double a25;
  double b25;
};

class zodiacalLight{

 private:
  vector<double> loglambda;
  vector<double> logI;
  felong_const fec;

  rbf wavelength;
  rbf position;

  //integrates using the trapezium rule
  double integrate(vector<double> x, vector<double> y);

  inline double felong(double elong)
  {
    return fec.a02 + fec.b02 * (elong>90.0?0.9:0.9+0.3*(elong-90)/60.0)
      + fec.a25 + fec.b25 * (elong>90.0?0.6:0.6+0.2*(elong-90)/60.0);
  };

 public:

  void outputSpectrum();

  zodiacalLight();
  ~zodiacalLight(){};

  void set_bandpass(vector<double> lambda, vector<double> eff);

  //counts per second over the bandpass
  inline double get_counts(double elong, double beta)
  {
    //get into the right range
    elong = fmod(abs(elong),360);
    elong = (elong<=180?elong:180-fmod(elong,180));
    beta = fmod(abs(beta),180);
    beta = (beta<=90?beta:90-fmod(beta,90));

    double eb[2]={elong,beta};
    
    return felong(elong)*exp(position.evaluate(eb,2));
  }

};

#define ZODIACAL_LIGHT_HEADER
#endif
