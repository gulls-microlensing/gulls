// VBMicrolensing v4.2 (2025)
//
// This code has been developed by Valerio Bozza (University of Salerno) and collaborators.
// Check the repository at https://github.com/valboz/VBMicrolensing
// for the newest version.
// Any use of the code for scientific publications should be acknowledged by a citation
// to the appropriate publication, as detailed in the repository page.
//
// The code relies on the root solving algorithm by Jan Skworon and Andy Gould
// described in Skowron & Gould arXiv:1203.1034.
// Please also cite this paper if specifically relevant in your scientific publication.
// The original Fortran code available on http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
// has been translated to C++ by Tyler M. Heintz and Ava R. Hoag (2017)
//
// GNU Lesser General Public License applies to all parts of this code.
// Please read the separate LICENSE.txt file for more details.

#pragma message("Including VBMicrolensingLibrary.h")
#ifndef __multilens
#define __multilens
#define __unmanaged

#define _L1 x1-((x1+a/2.0)/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*(x1-a/2.0)/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q) // Used in PlotCrits
#define _L2 x2-(x2/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*x2/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q)
#define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test
#define _J1c coefs[21]/((zc-coefs[20])*(zc-coefs[20]))+coefs[22]/(zc*zc) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J2 -2.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z))
#define _J3 6.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z*z))
#define _skew(p1,p2,q1,q2) p1*q2-p2*q1
#define _NP 200.0
#define __rsize 151
#define __zsize 101

#define _sign(x) ((x>0)? +1 : -1)

#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

class _curve;
class _sols;
class _theta;
class VBMcomplex;
struct annulus;


class VBMcomplex {
public:
	double re;
	double im;
	VBMcomplex(double a, double b) {re = a; im = b; }
	VBMcomplex(double a)           {re = a; im = 0; }
	VBMcomplex(void)               {re = 0; im = 0; }
};



class VBMicrolensing
{
	const int maxiter = 10000;
	const double pseudorandom[12] = { 0.5, 0.4, 0.3, 0.7, 0.6, 0.5, 0.4, 0.7, 0.5, 0.6, 0.3,
		0.3 };

	int *ndatasat;
	double **tsat,***possat;
	double Mag0, corrquad, corrquad2, safedist;
	double *dist_mp, *q;
	int nim0,n,n2,nnm1,nroots, nrootsmp, *nrootsmp_mp;
	VBMcomplex *zr, *zcr,**pmza, **pyaza , **ppmy , *pza , *pza2, **pmza2, *pdum , *ppy, *a, *s_offset, *pert,y,yc,*s;
	VBMcomplex *y_mp, *** pmza_mp, ** pza_mp, ***pyaza_mp, ***ppmy_mp, **ppy_mp, **zr_mp;
	VBMcomplex *zaltc, *J1, *J1c,**za,**za2; 
	VBMcomplex *coefs, **coefs_mp;
	VBMcomplex** a_mp, *s_sort;

	double *prodevs, *errs,err,L0f,Jacf;
	VBMcomplex *devs, *init, *centralimages, *newseeds,*grads,zf,S2f,*S2s,*S3s,*S4s;
	int lencentralimages,lennewseeds,ngoodold,ngood,iter,iter2;
	////
	double* good, * Jacs, rho, rho2, * m;
	double** m_mp, *q_sort;
	int *worst;
	double e,phi,phip,phi0,Om,inc,t0,d3,v3,GM,flagits;
	double Obj[3],rad[3],tang[3],t0old;
	double Eq2000[3],Quad2000[3],North2000[3]; 
	double ESPLout[__rsize][__zsize], ESPLin[__rsize][__zsize], ESPLoutastro[__rsize][__zsize], ESPLinastro[__rsize][__zsize];
	bool ESPLoff, multidark;
	double* LDtab, * rCLDtab, * CLDtab;
	double scr2, sscr2;
	int npLD;
	annulus *annlist;
	_curve **cprec, **cpres, **cfoll;
	double **A;
	
	void ComputeParallax(double, double, double *);
	double LDprofile(double r);
	double rCLDprofile(double tc, annulus*, annulus*);
	void initroot();
	int froot(VBMcomplex);
	bool checkroot(_theta *);

	void SetLensGeometry_spnp(int n, double* q, VBMcomplex* s);
	void SetLensGeometry_multipoly(int n, double* q, VBMcomplex* s);
	void initrootpoly();
	_curve *NewImages(VBMcomplex,VBMcomplex  *,_theta *);
	_curve *NewImages(_theta *);
	_curve *NewImagespoly(_theta *);
	_curve* NewImagesmultipoly(_theta*);
	double BinaryMagSafe(double s, double q, double y1, double y2, double rho, _sols** images);
	void OrderImages(_sols *,_curve *);
	void OrderMultipleImages(_sols *, _curve *);
	void cmplx_laguerre(VBMcomplex *, int, VBMcomplex *, int &, bool &);
	void cmplx_newton_spec(VBMcomplex *, int, VBMcomplex *, int &, bool &);
	void cmplx_laguerre2newton(VBMcomplex *, int, VBMcomplex *, int &, bool &, int);
	void solve_quadratic_eq(VBMcomplex &, VBMcomplex &, VBMcomplex *);
	void solve_cubic_eq(VBMcomplex &, VBMcomplex &, VBMcomplex &, VBMcomplex *);
	void polyproduct(VBMcomplex *p1, int n1, VBMcomplex *p2, int n2, VBMcomplex *pdest);
	void copypol(VBMcomplex *p1, int n1, VBMcomplex *pdest); 
	void change_n(int nn);
	void change_n_mp(int nn);
	void polycoefficients();
	void polycoefficients_multipoly();
	void polycritcoefficients(VBMcomplex eiphi);

public: 

//	bool testnewcoefs;

	void SetLensGeometry(int n, double* q, VBMcomplex *s);
	void SetLensGeometry(int n, double* pr);
	double MultiMag0(VBMcomplex y, _sols** Images);
	double MultiMag0(VBMcomplex y);
	double MultiMag0(double y1, double y2);
	double MultiMag(VBMcomplex y, double rho, double accuracy, _sols **Images);
	double MultiMag(VBMcomplex y, double rho, double accuracy);
	double MultiMag(VBMcomplex y, double rho);
	double MultiMag(double y1, double y2, double rho);
	double rootaccuracy;
	double samplingfactor;
	bool squarecheck;
	bool astrometry;

	static char ESPLtablefile[1024];
	static void SetESPLtablefile(char *instring) { strcpy(ESPLtablefile, instring); }
	double Tol,RelTol,a1,a2,t0_par;
	double mass_radius_exponent, mass_luminosity_exponent;
	int satellite,parallaxsystem,t0_par_fixed,nsat;
	int minannuli,maxannuli,nannuli,NPS,NPcrit;
	int newtonstep;
	double y_1,y_2,av, therr, astrox1,astrox2;
	double (*CumulativeFunction)(double r,double *LDpars);

// Critical curves and caustics calculation
	_sols* PlotCrit();
	_sols *PlotCrit(double a,double q);
// Initialization for parallax calculation
	void SetObjectCoordinates(char *Coordinates_file, char *Directory_for_satellite_tables);
	void SetObjectCoordinates(char *CoordinateString);
// Skowron & Gould root calculation
	void cmplx_roots_gen(VBMcomplex *, VBMcomplex *, int, bool, bool);
	void cmplx_roots_multigen(VBMcomplex*, VBMcomplex**, int, bool, bool);
// Bozza optimization
	int findimagepoly(int iroot);
	int findimagemultipoly(int iroot);

// Magnification calculation functions.

	double BinaryMag0(double s,double q,double y1,double y2, _sols **Images);
	double BinaryMag0(double s, double q, double y1, double y2);
	double BinaryMag(double s,double q,double y1,double y2,double rho,double accuracy, _sols **Images);
	double BinaryMag(double s,double q ,double y1,double y2,double rho,double accuracy);
	double BinaryMag2(double s, double q, double y1, double y2, double rho);
	double BinaryMagDark(double s, double q, double y1, double y2, double rho, double accuracy);
	void BinaryMagMultiDark(double s, double q, double y1, double y2, double rho, double *a1_list, int n_filters, double *mag_list, double accuracy);

// Limb Darkening control
	enum LDprofiles { LDlinear, LDquadratic, LDsquareroot, LDlog, LDuser };
	void SetLDprofile(double(*UserLDprofile)(double), int tablesampling);
	void SetLDprofile(LDprofiles);

// Method control
	enum class Method { Singlepoly, Multipoly, Nopoly};
	void SetMethod(Method);
        
//ESPL functions
	void LoadESPLTable(char *tablefilename);
	double ESPLMag(double u, double rho);
	double ESPLMag2(double u, double rho);
	double ESPLMagDark(double u, double rho);
	double PSPLMag(double u);


// New (v2) light curve functions, operating on arrays

	void PSPLLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void PSPLLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void ESPLLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void ESPLLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);

	void BinaryLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void BinaryLightCurveW(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void BinaryLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void BinaryLightCurveOrbital(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *sep_array, int np);
	void BinaryLightCurveKepler(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);

	void BinSourceLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void BinSourceLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void BinSourceLightCurveXallarap(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, double *sep_array, int np);
	void BinSourceExtLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceSingleLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceBinLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);

	void TripleLightCurve(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
	void TripleLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void LightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np, int nl);

// Old (v1) light curve functions, for a single calculation
	double PSPLLightCurve(double *parameters, double t);
	double PSPLLightCurveParallax(double *parameters, double t);
	double ESPLLightCurve(double *parameters, double t);
	double ESPLLightCurveParallax(double *parameters, double t);

	double BinaryLightCurve(double *parameters,double t);
	double BinaryLightCurveW(double *parameters, double t);
	double BinaryLightCurveParallax(double *parameters, double t);
	double BinaryLightCurveOrbital(double *parameters, double t);
	double BinaryLightCurveKepler(double* parameters, double t);

	double BinSourceLightCurve(double *parameters, double t);
	double BinSourceLightCurveParallax(double *parameters, double t);
	double BinSourceLightCurveXallarap(double *parameters, double t);
	double BinSourceExtLightCurve(double* parameters, double t);
	double BinSourceBinLensXallarap(double* parameters, double t);
	double BinSourceSingleLensXallarap(double* parameters, double t);
	double BinSourceBinLensPOX(double* parameters, double t);

	double TripleLightCurve(double *parameters, double t);
	double TripleLightCurveParallax(double* parameters, double t);
	// Constructor and destructor

	VBMicrolensing();
	~VBMicrolensing();

	private:
		LDprofiles curLDprofile;
		Method SelectedMethod;
};

//std::string VBMicrolensing::ESPLtablefile = "";

double VBDefaultCumulativeFunction(double r, double *a1);

struct annulus{
	double bin;
	double cum;
	double Mag;
	double err;
	double f;
	int nim;
	double LDastrox1, LDastrox2;
	annulus *prev,*next;
};



class _theta{
public: 
	double th,maxerr,Mag,errworst,astrox1,astrox2;
	int imlength;
	_theta *prev,*next;

	_theta(double);

};

class _thetas{
public: 
	_theta *first,*last;
	int length;

	_thetas(void);
	~_thetas(void);
	_theta *insert(double);
	void remove(_theta*);

};

class _point{
public:
	double x1;
	double x2;
	double parab,ds,dJ,Mag,err,parabastrox1, parabastrox2;;
	VBMcomplex d,J2;
	_theta *theta;
	_point(double ,double,_theta *);
	_point *next,*prev;
	double operator-(_point);
};

class _curve{
public:
	int length;
	_point *first,*last;
	_curve *next,*prev;
	_curve *partneratstart,*partneratend;
	double parabstart,Magstart,errstart, parabastrox1, parabastrox2;

	_curve(_point *);
	_curve(void);
	~_curve(void);

	_curve *divide(_point *);
	void drop(_point *);
	void append(double,double);
	void append(_point *);
	void prepend(double,double);
//	void prepend(_point *);
	_curve *join(_curve *);
	_curve *joinbefore(_curve *);
	_curve *reverse(void);
	double closest(_point *,_point **);
	double closest2(_point *,_point **);
	void complement(_point **,int,_point **,int);
};

class _sols{
public:
	int length;
	_curve *first,*last;

	_sols(void);
	~_sols(void);
	void drop(_curve *);
	void append(_curve *);
	void prepend(_curve *);
	void join(_sols *);	
};


#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 101

//////////////////////////////
//////////////////////////////
////////VBMcomplex methods and operators
//////////////////////////////
//////////////////////////////

inline double abs2(VBMcomplex z) {
	return (z.re * z.re + z.im * z.im);
}

inline double abs(VBMcomplex z) {
	return sqrt(z.re * z.re + z.im * z.im);
}

inline VBMcomplex conj(VBMcomplex z) {
	return VBMcomplex(z.re, -z.im);
}

inline VBMcomplex sqrt(VBMcomplex z) {
	double md = sqrt(z.re * z.re + z.im * z.im);
	return (md > 0) ? VBMcomplex((sqrt((md + z.re) / 2) * ((z.im > 0) ? 1 : -1)), sqrt((md - z.re) / 2)) : 0.0;
}

inline double real(VBMcomplex z) {
	return z.re;
}

inline double imag(VBMcomplex z) {
	return z.im;
}

inline VBMcomplex operator+(VBMcomplex p1, VBMcomplex p2) {
	return VBMcomplex(p1.re + p2.re, p1.im + p2.im);
}

inline VBMcomplex operator-(VBMcomplex p1, VBMcomplex p2) {
	return VBMcomplex(p1.re - p2.re, p1.im - p2.im);
}

inline VBMcomplex operator*(VBMcomplex p1, VBMcomplex p2) {
	return VBMcomplex(p1.re * p2.re - p1.im * p2.im, p1.re * p2.im + p1.im * p2.re);
}

inline VBMcomplex operator/(VBMcomplex p1, VBMcomplex p2) {
	double md = p2.re * p2.re + p2.im * p2.im;
	return VBMcomplex((p1.re * p2.re + p1.im * p2.im) / md, (p1.im * p2.re - p1.re * p2.im) / md);
}

inline VBMcomplex operator+(VBMcomplex z, double a) {
	return VBMcomplex(z.re + a, z.im);
}

inline VBMcomplex operator-(VBMcomplex z, double a) {
	return VBMcomplex(z.re - a, z.im);
}

inline VBMcomplex operator*(VBMcomplex z, double a) {
	return VBMcomplex(z.re * a, z.im * a);
}

inline VBMcomplex operator/(VBMcomplex z, double a) {
	return VBMcomplex(z.re / a, z.im / a);
}

inline VBMcomplex operator+(double a, VBMcomplex z) {
	return VBMcomplex(z.re + a, z.im);
}

inline VBMcomplex operator-(double a, VBMcomplex z) {
	return VBMcomplex(a - z.re, -z.im);
}

inline VBMcomplex operator*(double a, VBMcomplex z) {
	return VBMcomplex(a * z.re, a * z.im);
}

inline VBMcomplex operator/(double a, VBMcomplex z) {
	double md = z.re * z.re + z.im * z.im;
	return VBMcomplex(a * z.re / md, -a * z.im / md);
}

inline VBMcomplex operator+(VBMcomplex z, int a) {
	return VBMcomplex(z.re + a, z.im);
}

inline VBMcomplex operator-(VBMcomplex z, int a) {
	return VBMcomplex(z.re - a, z.im);
}

inline VBMcomplex operator*(VBMcomplex z, int a) {
	return VBMcomplex(z.re * a, z.im * a);
}

inline VBMcomplex operator/(VBMcomplex z, int a) {
	return VBMcomplex(z.re / a, z.im / a);
}

inline VBMcomplex operator+(int a, VBMcomplex z) {
	return VBMcomplex(z.re + a, z.im);
}

inline VBMcomplex operator-(int a, VBMcomplex z) {
	return VBMcomplex(a - z.re, -z.im);
}

inline VBMcomplex operator*(int a, VBMcomplex z) {
	return VBMcomplex(a * z.re, a * z.im);
}

inline VBMcomplex operator/(int a, VBMcomplex z) {
	double md = z.re * z.re + z.im * z.im;
	return VBMcomplex(a * z.re / md, -a * z.im / md);
}

inline VBMcomplex operator-(VBMcomplex z) {
	return VBMcomplex(-z.re, -z.im);
}

inline bool operator==(VBMcomplex p1, VBMcomplex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

inline bool operator!=(VBMcomplex p1, VBMcomplex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

inline VBMcomplex expcmplx(VBMcomplex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return VBMcomplex(r * cos(theta), r * sin(theta));
}

inline VBMcomplex cbrt(VBMcomplex z) {
	VBMcomplex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 0.333333333333);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	return 	VBMcomplex(r_cube * cos(theta_cube), r_cube * sin(theta_cube));
}

#endif
