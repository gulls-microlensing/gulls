#ifndef PARALLAX_HEADER

#include "ephem.h"
#include "coords.h"

using namespace std;

class parallax{

  //Tools for easily computing microlensing parallax. Relies on the 
  //orbits computed using ephem.

  //Conventions:
  //
  //  Units: AU, day, radians
  //  murel is the proper motion of the lens relative to the source
  //  _piE_ = (piEN, piEE) = (cos phipi, sin phipi) piE
  //
  //  When quantities differ between reference frames, the following
  //  subscripts will identify them, e.g. mu1_hg is the first 
  //  component of heliocentric proper motion in Galactic coordinates
  //  i.e., mu_l_helio
  //
  //    h = heliocentric (technically barycentric)
  //    r = reference-frame-centric
  //    o = observer-centric
  //    g = galactic
  //    c = celestial/equatorial
  //    e = ecliptic

  //

 private:

 public:
  
  int debug;
  static const int POSITION=1;  static const int LPOSITION=0;
  static const int REFFRAME=2;  static const int LREFFRAME=1;
  static const int ORBIT=4;     static const int LORBIT=2;
  static const int PROVIDED=8;  static const int LPROVIDED=3;
  static const int EPOCHS=16;   static const int LEPOCHS=4;
  static const int NESHIFTS=32; static const int LNESHIFTS=5;
  static const int TUSHIFTS=64; static const int LPARALLAX=6;
  static const int PLLXINITIALIZED=32768;   static const int LINITIALIZED=15;
  static const int allready=POSITION+ORBIT+REFFRAME+PROVIDED+EPOCHS
	+NESHIFTS+TUSHIFTS;
  static constexpr double daysinyr = 365.25;
  static constexpr double auday_to_kms = 1731.456837;

  //Coordinate transformations object
  coords c;

  //Reference frame. Flag: 2^0
  double tref; //Defined by an epoch and an orbit
  vector<orbitalElements>* oref; //orbit of the reference frame

  //      derived quantities
  vector<double> xref, vref, aref; //position, velocity and acceleration
                                   //of inertial reference frame origin at tref

  //Orbit. Flag: 2^1
  vector<orbitalElements>* orbit;

  //Position: Flag: 2^2
  double a, d, l, b, lam, bet;

  //Relative proper motion of the lens in the heliocentric frame
  //unit vectors. Flag: 2^3
  double mua_h, mud_h, mul_h, mub_h, mulam_h, mubet_h, murel_h;
  //in the reference frame
  double mua_r, mud_r, mul_r, mub_r, mulam_r, mubet_r, murel_r;

  //Projected velocities
  double vtilde_h, vtilde_r, v_ref;
  double vtilde_N_h, vtilde_E_h, vtilde_N_r, vtilde_E_r, v_ref_N, v_ref_E;
  

  //Unit vectors of parallax and components. Flag: 2^4
  double piEN, piEE; //Components of the parallax in the direction of lens motion in the reference frame
  //(axes are the same as for the ecliptic plane, just different origin)
  double piEll, piErp; //components of lens motion // and |_ to the Sun's acceleration vector in the reference frame
  double phi_llN; //angle between sun's acceleration vector (i.e. -a_ref) and N 
  double phi_pi; //angle between lens motion and N as measured in the reference frame
  double piE; //Magnitude
  double thetaE; //mas
  double rEtilde; //AU

  //uL parameters (defined at tref). Flag: 2^5
  double tE_h;
  double tE_r;
  //double u0_o, t0_o, tE_o;

  //Epochs. Flag: 2^6
  vector<double> epochs;

  //Status flags
  int status; //Do we have all the info we need?
  vector<int> progress; //for keeping track of what still needs 
                        //doing

  //vector<vector<double> > NEshift; //in AU
  vector<double> Nshift; //in AU
  vector<double> Eshift;
  vector<double> tshift; //in thetaE
  vector<double> ushift; //in thetaE
  vector<vector<double> > sslocation; //in AU in ecliptic Barycentric cartesian coordinates

  //constructors

  //Construct the class - reconstruct for a new orbit or viewing 
  //direction
  //eqecl allows it to work in the frame where NE is ecliptic (1) or equatorial(0)
  parallax();


  //Set-up functions

  void testinit(string fn, int status, int target);
  //int set_reference();
  int set_orbit(vector<orbitalElements>* orbit_);
  int set_radec(double ra, double dec, int deghr=0);
  int set_radec(vector<double> ra, vector<double> dec); //sexagesimal
  int set_lb(double l_, double b_, int deg=0); //deg=0 means in degrees
  //int set_pm_lb(double mul_, double mub_);
  //int set_pm_ad(double mua_, double mud_);
  //int set_piEpp(double piEll_, double piErp_);
  //int set_piE(double piE_);
  //int set_piENE(double piEN_, double piEE_);
  //int set_tE_h(double tE);
  //int set_tE_r(double tE);
  int load_epochs(double epochs_[], int n);
  int load_epochs(vector<double>* epochs_);

  int provide_murel_h_ad(double mua_h_, double mud_h_, double piE_, double thetaE_);
  int provide_murel_h_lb(double mul_h_, double mub_h_, double piE_, double thetaE_);
  int provide_observables_NE(double piEN_, double piEE_, double tE_r);
  int provide_observables_llrp(double piEll_, double piErp_, double tE_r);

  //void initialize();
  //void fit_reinit();
  void print_uninit();
  int compute_NEshifts();
  int compute_tushifts();
  int compute_tushifts(vector<int>* indices);
  int setup_reference_frame(double tref_, vector<orbitalElements>* oref_); //private
  //void compute_directions();


  void reset();

  /*  double tshift(int i)
  {
    return tshift[i];
  };
  double ushift(int i)
  {
    return ushift[i];
    };*/

  //  set_murel_hg(double mul, double mub);


  //For each epoch, compute the projection of the observer's orbit
  //on the plane perpendicular to the Sun-event direction
  //compute only once per orbit, direction and time-series
  //computeProjections();

  //Adjust the u (|_) and tau (//) coordinates of the source due to 
  //the parallax effect (depends on the piEE and piEN components 
  //that will need to be precomputed and stored at the class level
  //void utau_shift(int ep, double uin, double tauin, double* uout, 
  //		  double* tauout);
  //NOTE THAT THE DEFINITION OF // AND |_ DIFFER COMPLETELY FROM THE
  //LITERATURE DEFINITION, WHERE THEY CORRESPOND TO // AND |_ TO THE
  //SOLAR ACCELERATION VECTOR. ll AND rp INDICATE THE SUN ACCELERATION
  //DIRECION VALUES

  //Compute the timescale and relative proper motion of the event
  //that will be observed by the reference observer (given the 
  //heliocentric values), as well as the EN->//|_ transformation
  //matrix
  //  observed_tEmurel();
  

  //Compute the position of the observer at each epoch in the 
  //heliocentric frame where n points towards the ra,dec, and E is 
  //eastward, N northward in AU epochs should be in HJD
  //void computePositions(double ra, double dec, orbitalElements* ephem vector<double>* epochs, vector<double>* E, vector<double>* N, vector<double>* n);

};

#define PARALLAX_HEADER
#endif
