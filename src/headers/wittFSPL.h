#ifndef WITTPLFS

#include<vector>

#include"cd.h"
#include"constants.h"
#include"src_cld.h"


using namespace std;


////////////////////////////////////////////////////////////
//
//    Finite source point lens magnification
//      taken from Witt & Mao, 1994, ApJ, 430, 505
//      http://adsabs.harvard.edu/abs/1994ApJ...430..505W
//
///////////////////////////////////////////////////////////

struct wfspl_parameters {double zs; double rs; double Kx; double Ky;};
struct innerParams {src_cld* src; double zs; double costheta;}; 
struct outerParams {src_cld* src; double zs;}; 

static const double gphOffX[13]={0.0, 0.5,0.0,-0.5,0.0, 1.0,0.0,-1.0,0.0, oosqrt2,-oosqrt2,-oosqrt2,oosqrt2};

static const double gphOffY[13]={0.0, 0.0,0.5,0.0,-0.5, 0.0,1.0,0.0,-1.0, oosqrt2,oosqrt2,-oosqrt2,-oosqrt2};

double wittFSMagnification(double zsc, double rs);

inline double wittFSMagnification(cd zsc_, double rs)
{
  return wittFSMagnification(abs(zsc_),rs);
};

void dmududr_witt(double u, double rs, double* dmudu, double* dmudr, double delta=1e-6);

inline double wfsmag(cd zsc_, double rs)
{
  return wittFSMagnification(zsc_, rs);
};

inline double wfsmag(double zsc_, double rs)
{
  return wittFSMagnification(zsc_, rs);
};

double fsplMagnification(double zs, src_cld* src);
double muFunc(double phi, void *parameter);
double gphmag(double zs, src_cld* src, double& Asecond);

double leeFSPLmagnification(double zs, src_cld* src);
double leeOuter(double x, void * params);
double leeInner(double x, void * params);

double leeFSMagnification(double u, double r);
double lee_f(double theta, double u, double r);
double lee_fg(double theta, double u, double r, double asinpu);

#define WITTPLFS
#endif /* WITTPLFS */
