//        Functions for estimating the intensity of the zodiacal light
//                        in a user defined bandpass
//
//        Based on the tabulations of Leinert et al A&ASS 127,1 (1998)
//
//                     Copyright (C) 2012 Matthew Penny
//                      penny@astronomy.ohio-state.edu

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

  //void set_bandpass(vector<double> lambda, vector<double> eff);
  //nu should be in PHz
  void set_bandpass(vector<double> nu, vector<double> eff);
  void set_bandpass(string filename);

  inline double get_flux(double elong, double beta)
  {
    //get into the right range
    elong = fmod(abs(elong),360);
    elong = (elong<=180?elong:180-fmod(elong,180));
    beta = fmod(abs(beta),180);
    beta = (beta<=90?beta:90-fmod(beta,90));

    double eb[2]={elong,beta};
    //cout << "felong: " << felong(elong) << " " << exp(position.evaluate(eb,2)) << " " << position.evaluate(eb,2) << " " << elong << " " << beta << endl;
    
    return felong(elong)*exp(position.evaluate(eb,2));
  }

  inline double get_mag(double elong, double beta)
  {
    return 8.9-2.5*log10(get_flux(elong,beta));
  }
};

#define ZODIACAL_LIGHT_HEADER
#endif
