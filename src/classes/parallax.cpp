#include<cmath>
#include<vector>
#include<iostream>

#include "parallax.h"
#include "ephem.h"
#include "coords.h"
#include "constants.h"
#include "integerPowers.h"

parallax::parallax()
{
  status=0;
  debug=0;
  progress = vector<int>(16,0);
}

void parallax::testinit(string fn, int status, int target)
{
  if(debug) cout << __FUNCTION__ << endl;
  
  if((status&(target-1)) != (target-1))
    {
      cerr << __FILE__ << ": " << fn << ": Error: Some element of the parallax computation is uninitialized (all those unininitialized before this function must be called prior to running this:" << endl;
      print_uninit();
      exit(1);
    }

}


/* obsolete
void parallax::initialize()
{
  setup_reference_frame();
  compute_NEshifts();
  fit_reinit();
}
*/
/*
void parallax::fit_reinit()
{
  //Reinitialize assuming only tE and parallax vector has changed
  compute_directions();
  status |= PLLXINITIALIZED;
  compute_tushifts();
}
*/




void parallax::reset()
{
  if(debug) cout << __FUNCTION__ << endl;
  epochs.clear();
  //NEshift.clear();
  Nshift.clear();
  Eshift.clear();
  tshift.clear();
  ushift.clear();
  progress = vector<int>(16,0);
  status = 0;
}





void parallax::print_uninit()
{
  if(debug) cout << __FUNCTION__ << endl;
  if(!(status & POSITION)) cerr << "POSITION Uninitialized" << endl;
  if(!(status & REFFRAME)) cerr << "REFFRAME Uninitialized" << endl;
  if(!(status & ORBIT)) cerr << "ORBIT Uninitialized" << endl;
  if(!(status & PROVIDED)) cerr << "PROVIDED Uninitialized - run a provide... function" << endl;
  if(!(status & EPOCHS)) cerr << "EPOCHS Uninitialized" << endl;
  if(!(status & NESHIFTS)) cerr << "NESHIFTS Uninitialized" << endl;
  if(!(status & TUSHIFTS)) cerr << "TUSHIFTS Uninitialized" << endl;
}





int parallax::compute_NEshifts()
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,NESHIFTS);
  
  vector<double> x(3), xp;
  //NEshift = vector<vector<double> >(epochs.size(),vector<double>(2,0));
  Nshift = vector<double>(epochs.size());
  Eshift = vector<double>(epochs.size());
  tshift = vector<double>(epochs.size());
  ushift = vector<double>(epochs.size());
  sslocation = vector<vector<double> >(epochs.size(),vector<double>(3,0));

  for(int i=0;i<int(epochs.size());i++)
    {
      //compute the shift relative to the reference frame
      //x = vector<double>(3,0);
      x[0]=0; x[1]=0; x[2]=0;
      //orbit->viewfrom(epochs[i],a,d,&x);
      for(int j=0;j<int(orbit->size());j++)
	{
	  (*orbit)[j].viewfrom(epochs[i],a,d,&xp);
	  for(int k=0;k<int(x.size());k++)
	    {
	      x[k] += xp[k];
		}
	  sslocation[i][0] += (*orbit)[j].xecl;
	  sslocation[i][1] += (*orbit)[j].yecl;
	  sslocation[i][2] += (*orbit)[j].zecl;
	    

	  //figure out how to print out the x,y,z coordinates in the ecliptic plane here
	}
      //NEshift[i][0] = x[0] - xref[0] - (epochs[i]-tref)*vref[0];
      //NEshift[i][1] = x[1] - xref[1] - (epochs[i]-tref)*vref[1];
      Nshift[i] = x[0] - xref[0] - (epochs[i]-tref)*vref[0];
      Eshift[i] = x[1] - xref[1] - (epochs[i]-tref)*vref[1];
    }
  status |= NESHIFTS;
  return status;
}





int parallax::compute_tushifts()
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,TUSHIFTS);
 
  //tushift = vector<vector<double> >(epochs.size(),vector<double>(2,0));

  double cs= cos(phi_pi); 
  double sn= sin(phi_pi);

  //cout << "epochs.size() " << epochs.size() << endl;

  for(int i=0;i<int(epochs.size());i++)
    {
      //Convert the shift in the observer plane to a shift in the 
      //source position

      //tshift[i] = -piE * ( NEshift[i][0]*cs + NEshift[i][1]*sn);
      //ushift[i] = -piE * (-NEshift[i][0]*sn + NEshift[i][1]*cs);
      tshift[i] = -piE * ( Nshift[i]*cs + Eshift[i]*sn);
      ushift[i] = -piE * (-Nshift[i]*sn + Eshift[i]*cs);
    }
  status |= TUSHIFTS;
  return status;
}


/*
int parallax::compute_tushifts(vector<vector<double> >* tushift)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,TUSHIFTS);
  

  if(tushift->size()!=epochs.size())
    {
      if(tushift->size()==0)
	{
	  tushift.resize(epochs.size(),vector<double>(2,0));
	}
      else
	{
	  cerr << "compute_tushifts passed pointer to vector of the wrong size, needs to be " << epochs.size() << " received " << tushifts->size() << endl;
	  exit(1);
	}
    }

  double cs= cos(phi_pi); 
  double sn= sin(phi_pi);

  for(int i=0;i<int(epochs.size());i++)
    {
      //Convert the shift in the observer plane to a shift in the 
      //source position

      (*tushift)[i][0] = -piE * ( NEshift[i][0]*cs + NEshift[i][1]*sn);
      (*tushift)[i][1] = -piE * (-NEshift[i][0]*sn + NEshift[i][1]*cs);
    }
  status |= TUSHIFTS;
  return status;
}
*/



int parallax::setup_reference_frame(double tref_, vector<orbitalElements>* oref_) //private
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,REFFRAME);

  tref=tref_;
  oref = oref_;

  //Set up the reference frame
  //oref->viewfrom(tref,a,d,&xref);
  //oref->velocity(tref,a,d,&vref);
  //oref->acceleration(tref,a,d,&aref);
  
  //reset the orbits
  xref=vref=aref=vector<double>(3,0);

  //The components of these vectors are 0=N, 1=E, 2=normal

  //add on orbit and purturbations, and compute the position, velocity, and acceleration in AU/day
  vector<double> xp, vp, ap;
  for(int i=0;i<int(oref->size());i++)
    {
	  //if(debug) cout << "obs " << i << " viewfrom " << tref << " " << a << " " << d << endl;
      (*oref)[i].viewfrom(tref,a,d,&xp);
      (*oref)[i].velocity(tref,a,d,&vp);
      (*oref)[i].acceleration(tref,a,d,&ap);
      for(int j=0;j<int(xp.size());j++) xref[j]+=xp[j];
      for(int j=0;j<int(vp.size());j++) vref[j]+=vp[j];
      for(int j=0;j<int(ap.size());j++) aref[j]+=ap[j];
    }

  status |= REFFRAME;
  return status;
}

int parallax::provide_murel_h_ad(double mua_h_, double mud_h_, double piE_, double thetaE_)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,PROVIDED);

  piE = piE_;
  thetaE = thetaE_;
  
  //cout << "Compute the heliocentric proper motion unit vectors in each coordinate system" << endl;
  murel_h = qAdd(mua_h_, mud_h_);
  mua_h = mua_h_/murel_h;
  mud_h = mud_h_/murel_h;
  c.muad2lb(a,d,mua_h,mud_h,&mul_h,&mub_h);
  c.muad2ecl(a,d,mua_h,mud_h,&mulam_h,&mubet_h);
  tE_h = thetaE/murel_h * daysinyr;
  

  //cout << "Compute omega_E_helio and geo (omega_E = 1/t_E)" << endl;
  vector<double> wh(2), wr(2), wref(2);//, wo(2); //inverse timescale vector 
  //wr = (piEN,piEE)/tE is the lens velocity in Einstein radii per day (piEN, piEE is the proper motion unit vector
  //in the reference frame)
  //wref = apparent motion of the reference frame in the same units, as viewed 
  //from the sun

  rEtilde = 1.0/piE; //in AU

  //cout << "Velocity of the reference frame in Einstein ring radii per day in ecliptic coordinates in the heliocentric frame" << endl;
  wref[0] = vref[1]*piE;
  wref[1] = vref[0]*piE; //Yes these should be switched
  v_ref_E = vref[1] * auday_to_kms;
  v_ref_N = vref[0] * auday_to_kms;
  v_ref = qAdd(v_ref_N,v_ref_E);
  //cout << "gal wref: " << wref[0] << " " << wref[1] << endl;

  //cout << "Velocity of the lens in heliocentric frame in ecliptic coordinates (in Einstein radii per day)" << endl;
  wh[0] = mulam_h/tE_h;
  wh[1] = mubet_h/tE_h;
  vtilde_E_h = wh[0] * rEtilde * auday_to_kms;
  vtilde_N_h = wh[1] * rEtilde * auday_to_kms;
  vtilde_h = qAdd(vtilde_N_h,vtilde_E_h);
  //cout << "gal wh: " << wh[0] << " " << wh[1] << endl;

  //cout << "Velocity of the lens in the reference frame in ecliptic coordinates (in Einstein radii per day)" << endl;
  wr[0] = wh[0] - wref[0];
  wr[1] = wh[1] - wref[1];
  vtilde_E_r = wr[0] * rEtilde * auday_to_kms;
  vtilde_N_r = wr[1] * rEtilde * auday_to_kms;
  vtilde_r = qAdd(vtilde_N_r,vtilde_E_r);
  //cout << "gal wr: " << wr[0] << " " << wr[1] << endl;

  //cout << "Timescale in the reference frame (this is what's measured from the lightcurve)" << endl;
  tE_r = 1.0/qAdd(wr[0],wr[1]);
  murel_r = thetaE/tE_r * daysinyr;

  //cout << "unit vector in ecliptic coordinates is in the same direction as proper motion" << endl;
  mulam_r = wr[0]*tE_r;
  mubet_r = wr[1]*tE_r;

  //cout << "Convert to the other coordinates, still in the reference frame" << endl;
  c.muecl2ad(a,d,mulam_r,mubet_r,&mua_r,&mud_r);
  c.muad2lb(a,d,mua_r,mud_r,&mul_r,&mub_r);
  
  //cout << "For parallax components, the unit vector is the same, so just multiply by piE" << endl;
  piEN = mubet_r * piE;
  piEE = mulam_r * piE;
  phi_pi = atan2(piEE,piEN);
  
  //cout << "Calculate proper motion components in llrp (// and perpendicular to the Sun's acceleration vector) - these should be the principal components of the error ellipse for annual parallax" << endl;
  phi_llN = atan2(-aref[1],-aref[0]);
  //cout << "phi_llN " << phi_llN << endl;
  piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
  piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);
  status |= PROVIDED;
  return status;
}

int parallax::provide_murel_h_lb(double mul_h_, double mub_h_, double piE_, double thetaE_)
{
  if(debug) cout << __FUNCTION__ << endl;
  double mua_h_, mud_h_;
  c.mulb2ad(l,b,mul_h_,mub_h_,&mua_h_,&mud_h_);
  provide_murel_h_ad(mua_h_,mud_h_,piE_,thetaE_);
  return status;
}

int parallax::provide_observables_NE(double piEN_, double piEE_, double tE_r_)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,PROVIDED);
 
  //cout << "Compute the proper motion unit vectors in the reference frame in different coordinates" << endl;

  vector<double> wh(2), wr(2), wref(2);//, wo(2); //inverse timescale vector 

  piE = qAdd(piEN_,piEE_);
  piEN = piEN_; 
  piEE = piEE_;
  phi_pi = atan2(piEE,piEN);
  tE_r = tE_r_;

  rEtilde = 1.0/piE; //in AU

  //cout << "Calculate proper motion components in llrp (// and perpendicular to the Sun's acceleration vector) - these should be the principal components of the error ellipse for annual parallax" << endl;
  phi_llN = atan2(-aref[1],-aref[0]);
  //cout << "phi_llN " << phi_llN << endl;
  piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
  piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);

  //cout << "The unit vector is the same" << endl;
  mubet_r = piEN/piE;
  mulam_r = piEE/piE;

  //cout << "Convert to the other coordinates, still in the reference frame" << endl;
  c.muecl2ad(a,d,mulam_r,mubet_r,&mua_r,&mud_r);
  c.muad2lb(a,d,mua_r,mud_r,&mul_r,&mub_r);

  //cout << "We don't know thetaE, so let's just set mu_rel as if theta_E were 1.0 mas" << endl;
  murel_r = 1.0/tE_r * daysinyr;

  //cout << "Compute the vecocity of the observer in the reference frame in ecliptic coordinates (in Einstein radii per day)" << endl;
  wr[0] = mulam_r/tE_r;
  wr[1] = mubet_r/tE_r;
  vtilde_E_r = wr[0] * rEtilde * auday_to_kms;
  vtilde_N_r = wr[1] * rEtilde * auday_to_kms;
  vtilde_r = qAdd(vtilde_N_r,vtilde_E_r);
  //cout << "obs wr: " << wr[0] << " " << wr[1] << endl;

  //cout << "Velocity of the reference frame in Einstein ring radii per day in ecliptic coordinates in the heliocentric frame" << endl;
  wref[0] = vref[1]*piE; //yes they should be switched
  wref[1] = vref[0]*piE;
  v_ref_E = vref[1] * auday_to_kms;
  v_ref_N = vref[0] * auday_to_kms;
  v_ref = qAdd(v_ref_N,v_ref_E);
  //cout << "obs wref: " << wref[0] << " " << wref[1] << endl;

  //cout << "Velocity of the lens in the heliocentric frame in ecliptic coordinates (in Einstein radii per day)" << endl;
  wh[0] = wr[0] + wref[0];
  wh[1] = wr[1] + wref[1];
  vtilde_E_h = wh[0] * rEtilde * auday_to_kms;
  vtilde_N_h = wh[1] * rEtilde * auday_to_kms;
  vtilde_h = qAdd(vtilde_N_h,vtilde_E_h);
  //cout << "obs wh: " << wh[0] << " " << wh[1] << endl;

  //cout << "Unit vector for proper motion in the heliocentric frame" << endl;
  mulam_h = wh[0]/qAdd(wh[0],wh[1]);
  mubet_h = wh[1]/qAdd(wh[0],wh[1]);
  c.muecl2ad(a,d,mulam_h,mubet_h,&mua_h,&mud_h);
  c.muad2lb(a,d,mua_h,mud_h,&mul_h,&mub_h);
  
  //cout << "What's its magnitude" << endl;

  tE_h = 1.0/qAdd(wh[0],wh[1]); //This assumes thetaE is 1.0 mas
  murel_h = 1.0/tE_h * daysinyr; //This assumes thetaE is 1.0 mas
  status |= PROVIDED;
  return status;
}


int parallax::provide_observables_llrp(double piEll_, double piErp_, double tE_r_)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,PROVIDED);
  phi_llN = atan2(-aref[1],-aref[0]);
  double piEN_ = piEll*cos(phi_llN) - piErp*sin(phi_llN);
  double piEE_ = piEll*sin(phi_llN) + piErp*cos(phi_llN);

  provide_observables_NE(piEN_, piEE_, tE_r_);
  return status;
}

//This is obsolete and buggy (potentially)
/*
void parallax::compute_directions()
{
  //Add or subtract the positions and velocities of the various frames 
  //to get the event parameters in other frames.

  int exitnow=0;

  if(!(status&REFFRAME))
    {
      cerr << __FILE__ << ": Error: Reference frame not set" << endl;
      exitnow=1;
    }

  if(!(status&ORBIT))
    {
      cerr << __FILE__ << ": Error: Orbit not set (either through piE components or proper motions)." << endl;
      exitnow=1;
    }

  if(!(status&POSITION))
    {
      cerr << __FILE__ << ": Error: Event position not set" << endl;
      exitnow=1;
    }

  if(!(status&SOURCEDIR))
    {
      cerr << __FILE__ << ": Error: Source direction not set (either through piE components or proper motions)." << endl;
      exitnow=1;
    }

  if(!(status&PARAMS))
    {
      cerr << __FILE__ << ": Error: Event parameters not set." << endl;
      exitnow=1;
    }

  if(exitnow) exit(1);

  vector<double> x, v, xp, vp; //position and velocity of the object at tref

  //h = 1; r = 2;

  vector<double> wh(2), wr(2), wref(2);//, wo(2); //inverse timescale vector 
  //wr = (piEN,piEE)/tE is the lens velocity in Einstein radii per day (piEN, piEE is the proper motion unit vector
  //in the reference frame)
  //wref = apparent motion of the reference frame in the same units, as viewed 
  //from the sun

  wref[0] = vref[0]*piE;
  wref[1] = vref[1]*piE;
  //cout << "wref: " << wref[0] << " " << wref[1] << endl;

  //Calculate proper motion components in NE (ecliptic) and llrp (// and perpendicular to the acceleration vector)
  phi_llN = atan2(-aref[1],-aref[0]);
  //cout << "phi_llN " << phi_llN << endl;

  //set_pm_lb(mul_,mub_)      progress[LSOURCEDIR]=1;
  //set_pm_ad(mua_,mud_)      progress[LSOURCEDIR]=1;
  //set_piEpp(piEll_,piErp_)  progress[LSOURCEDIR]=2;
  //set_piENE(piEN_,piEE_)    progress[LSOURCEDIR]=3;

  //cerr << progress[LPARAMS] << " " << progress[LSOURCEDIR] << endl;

  int eqecl=1;

  if(progress[LPARAMS]==1 && progress[LSOURCEDIR]==1)
    {
	  if(eqecl==1)
		{
		  wh[0] = mulam_h/tE_h;
		  wh[1] = mubet_h/tE_h;
		}
	  else
		{
		  wh[0] = mud_h/tE_h; 
		  wh[1] = mua_h/tE_h;
		}

	  //cout << "wh: " << wh[0] << " " << wh[1] << endl;
      wr[0] = wh[0] - wref[0];
      wr[1] = wh[1] - wref[1];
	  //cout << "wr: " << wr[0] << " " << wr[1] << endl;
      tE_r = 1.0/qAdd(wr[0],wr[1]);
      piEN = wr[0]*tE_r;
      piEE = wr[1]*tE_r;
      piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
      piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);
	  //cout << "1 1: tE_r piEN piEE piEll piErp " << tE_r << " " << piEN << " " << piEE << " " << piEll << " " << piErp << endl;
      progress[LSOURCEDIR]=2;
    }
  else if(progress[LPARAMS]==2 && progress[LSOURCEDIR]>1)
    {
      if(progress[LSOURCEDIR]==2)
		{
		  piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
		  piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);
		  //cout << progress[LPARAMS] << " " << progress[LSOURCEDIR] << " piEll piErp: " << piEll << " " << piErp << endl;
		}
      else if(progress[LSOURCEDIR]==3)
		{
		  piEN = piEll*cos(phi_llN) - piErp*sin(phi_llN);
		  piEE = piEll*sin(phi_llN) + piErp*cos(phi_llN);
		  //cout << progress[LPARAMS] << " " << progress[LSOURCEDIR] << " piEll piErp: " << piEN << " " << piEE << endl;
		}
      wr[0] = piEN/tE_r;
      wr[1] = piEE/tE_r;
      wh[0] = wr[0] + wref[0];
      wh[1] = wr[1] + wref[1];
      tE_h = 1.0/qAdd(wh[0],wh[1]);
	  //cout << progress[LPARAMS] << " " << progress[LSOURCEDIR] << " wr wh tE_h: " << wr[0] << " " << wr[1] << " " << wh[0] << " " << wh[1] << " " << tE_h << endl;
	  if(eqecl==1)
		{
		  mulam_h = wh[0]*tE_h;
		  mubet_h = wh[0]*tE_h;
		  c.muecl2ad(a,d,mulam_h,mubet_h,&mua_h,&mud_h);
		  //cout << "AAA changed mua_h mud_h mulam_h mubet_h" << mua_h << " " << mub_h << " " << mulam_h << " " << mubet_h << endl;
		}
	  else
		{
		  mud_h = wh[0]*tE_h;
		  mua_h = wh[1]*tE_h;
		  c.muad2ecl(a,d,mua_h,mud_h,&mulam_h,&mubet_h);
		  //cout << "BBB changed mua_h mud_h mulam_h mubet_h" << mua_h << " " << mub_h << " " << mulam_h << " " << mubet_h << endl;
		}
      c.muad2lb(a,d,mua_h,mud_h,&mul_h,&mub_h);
	  //cout << progress[LPARAMS] << " " << progress[LSOURCEDIR] << " mua_h mud_h mulam_h mubet_h mul_h mub_h " << mua_h << " " << mud_h << " " << mulam_h << " " << mubet_h << " " << mul_h << " " << mub_h << endl; 
    }
  else
    {
      cerr << __FILE__ << ": " << string(__FUNCTION__) << ": Error: Must specify parameters and source directions as heliocentric or reference-frame-centric" << endl;
    }

  phi_pi = atan2(piEE,piEN);
  //cout << progress[LPARAMS] << " " << progress[LSOURCEDIR] << " phi_pi " << phi_pi << endl;
  //cout << "l b ra dec " << l << " " << b << " " << a << " " << d << endl;
  
}
*/

/* obsolete
int parallax::set_reference(double tref_, vector<orbitalElements>* oref_)
{
  testinit(string(__FUNCTION__),status,REFFRAME);
  tref=tref_;
  oref = oref_;

  status |= REFFRAME;
  return status;
}
*/

int parallax::set_orbit(vector<orbitalElements>* orbit_)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,ORBIT);
  orbit = orbit_;
  status |= ORBIT;
  return status;
}

int parallax::set_radec(double ra, double dec, int deghr)
{
  if(debug) cout << __FUNCTION__ << endl;
  //deghr = 0 - both radians
  //      = 1 - both degrees
  //      = 2 - hr and degrees

  if(deghr>2)
    {
      cerr << __FILE__ << ": "<< string(__FUNCTION__) << ": Error setting position - bad angle unit code (" << deghr << ")" << endl;
      exit(1);
    }

  if(deghr==0)
    {
      a = c.fold(ra,0,twoPi); 
      d = c.fold(dec,0,twoPi);
    }
  else if(deghr==1)
    {
      a = c.fold(c.deg2rad(ra),0,twoPi); 
      d = c.fold(c.deg2rad(dec),0,twoPi);
    }
  else
    {
      a = c.fold(c.hrs2rad(ra),0,twoPi); 
      d = c.fold(c.deg2rad(dec),0,twoPi);
    }

  c.ad2lb(a,d,&l,&b);
  c.ad2ecl(a,d,&lam,&bet);

  status = POSITION;
  return status;
}

int parallax::set_radec(vector<double> ra, vector<double> dec)
{
  if(debug) cout << __FUNCTION__ << endl;
  //deghr = 0 - both radians
  //      = 1 - both degrees
  //      = 2 - hr and degrees

  if(int(ra.size())!=3||int(dec.size())!=3)
    {
      cerr << __FILE__ << ": "<< string(__FUNCTION__) << ": Error setting position - input sexagesimal vectors (ra_hr,ra_min,ra_sec), (dec_deg,dec_min,dec_sec)" << endl;
      exit(1);
    }

  double a_,d_;
  int sgn_a, sgn_d;
  sgn_a = (ra[0]<0?-1:1);
  sgn_d = (dec[0]<0?-1:1);
  a_ = c.deg2rad((ra[0] + sgn_a*(ra[1]/60.0+ra[2]/3600.0))*15.0);
  d_ = c.deg2rad(dec[0] + sgn_d*(dec[1]/60.0+dec[2]/3600.0));

 
  a = c.fold(a_,0,twoPi); 
  d = c.fold(d_,0,twoPi);

  c.ad2lb(a,d,&l,&b);
  c.ad2ecl(a,d,&lam,&bet);

  status |= POSITION;
  return status;
}

int parallax::set_lb(double l_, double b_, int deg)
{
  //deg = 0 - both degrees
  //    = 1 - both radians

  if(deg>1)
    {
      cerr << __FILE__ << ": "<< string(__FUNCTION__) << ": Error setting position - bad angle unit code (" << deg << ")" << endl;
      exit(1);
    }

  if(deg==0)
    {
      l = c.fold(c.deg2rad(l_),0,twoPi); 
      b = c.fold(c.deg2rad(b_),0,twoPi);
    }
  else if(deg==1)
    {
      l = c.fold(l_,0,twoPi); 
      b = c.fold(b_,0,twoPi);
    }

  c.lb2ad(l,b,&a,&d);
  c.ad2ecl(a,d,&lam,&bet);
  
  status = POSITION;
  return status;
}

/* obsolete
int parallax::set_pm_lb(double mul_, double mub_)
{
  //set the heliocentric proper motion of the lens relative to the source in 
  //Galactic coordinates
  //assume mu_l = mu_l cos b ALWAYS
  
  if(!(status & POSITION))
    {
      cerr << __FILE__ << ": " << string(__FUNCTION__) << "Must set position in order to do proper motion unit conversions" << endl;
      exit(1);
    }
    
  double mu = qAdd(mul_,mub_);
  mul_h = mul_/mu;
  mub_h = mub_/mu;
  c.mulb2ad(l,b,mul_h,mub_h,&mua_h,&mud_h);
  c.muad2ecl(a,d,mua_h,mud_h,&mulam_h,&mubet_h);
  //cout << "helio mulb set, muad and muecl calc " << mua_h << " " << mud_h << " " << mulam_h << " " << mubet_h << endl;

  progress[LSOURCEDIR]=1;
  status |= SOURCEDIR;
  return status;
  
}

int parallax::set_pm_ad(double mua_, double mud_)
{
  //set the heliocentric proper motion of the lens relative to the source in 
  //equatorial coordinates
  //assume mu_a = mu_a cos d ALWAYS
  
  if(!(status & POSITION))
    {
      cerr << __FILE__ << ": " << string(__FUNCTION__) << "Must set position in order to do proper motion unit conversions" << endl;
      exit(1);
    }
    
  double mu = qAdd(mua_,mud_);
  mua_h = mua_/mu;
  mud_h = mud_/mu;
  c.muad2lb(a,d,mua_h,mud_h,&mul_h,&mub_h);
  c.muad2ecl(a,d,mua_h,mud_h,&mulam_h,&mubet_h);
  //cout << "helio muad set, mulb and muecl calc " << mua_h << " " << mud_h << " " << mulam_h << " " << mubet_h << endl;

  progress[LSOURCEDIR]=1;
  status |= SOURCEDIR;
  return status;
  
}
/*

/* obsolete
int parallax::set_piEpp(double piEll_, double piErp_)
{
  piE = qAdd(piEll_,piErp_);
  piEll = piEll_/piE;
  piErp = piErp_/piE;
  
  progress[LSOURCEDIR]=3;
  status |= PARALLAX;
  status |= SOURCEDIR;
  return status;
}
*/
 /* obsolete
int parallax::set_piE(double piE_)
{
  piE = piE_;

  status |= PARALLAX;
  return status;
}
 */
  /* obsolete
int parallax::set_piENE(double piEN_, double piEE_)
{
  piE = qAdd(piEN_,piEE_);
  piEN=piEN_/piE;
  piEE=piEE_/piE;
  
  progress[LSOURCEDIR]=2;
  status |= PARALLAX;
  status |= SOURCEDIR;
  return status;
}
  */
   /* obsolete
int parallax::set_tE_h(double tE)
{
  tE_h = tE;

  progress[LPARAMS]=1;
  status |= PARAMS;
  return status;
}
   */
   
	/* obsolete
int parallax::set_tE_r(double tE)
{
  tE_r = tE;

  progress[LPARAMS]=2;
  status |= PARAMS;
  return status;
}
	*/

/*int parallax::set_params_o(double t0, double tE, double u0)
{
  t0_o = t0;
  tE_o = tE;
  u0_o = u0;

  progress[LPARAMS]=3;
  status |= PARAMS;
  return status;
  }*/

int parallax::load_epochs(double epochs_[], int n)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,EPOCHS);
  
  epochs.resize(n);
  for(int i=0;i<n;i++)
    {
      epochs[i]=epochs_[i];
    }

  if(n>0) status |= EPOCHS;
  return status;
}

int parallax::load_epochs(vector<double>* epochs_)
{
  if(debug) cout << __FUNCTION__ << endl;
  testinit(string(__FUNCTION__),status,EPOCHS);
  int n=epochs_->size();
  epochs.resize(n);
  for(int i=0;i<n;i++)
    {
      epochs[i]=(*epochs_)[i];
    }

  if(n>0) status |= EPOCHS;
  return status;
}







/*void parallax::utau_shift(int ep, double uin, double tauin, double* uout, double* tauout)
{
  
  *tauout = tauin - piE*(EN2llrp00*obsshift[ep][0] 
			 + EN2llrp01*obsshift[ep][1]);
  *uout = uin - piE*(EN2llrp10*obsshift[ep][0] 
			 + EN2llrp11*obsshift[ep][1]);
			 }*/
