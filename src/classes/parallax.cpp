#include<cmath>
#include<vector>

#include "parallax.h"
#include "ephem.h"
#include "coords.h"
#include "constants.h"
#include "integerPowers.h"

parallax::parallax()
{
  status=0;
  progress = vector<int>(16,0);
}

void parallax::initialize()
{
  setup_reference_frame();
  compute_NEshifts();
  fit_reinit();
}

void parallax::fit_reinit()
{
  //Reinitialize assuming only tE and parallax vector has changed
  compute_directions();
  status |= PLLXINITIALIZED;
  compute_tushifts();
}

void parallax::reset()
{
  epochs.clear();
  NEshift.clear();
  tushift.clear();
  progress = vector<int>(16,0);
  status = 0;
}

void parallax::print_uninit()
{
  if(!(status & REFFRAME)) cerr << "REFFRAME Uninitialized" << endl;
  if(!(status & ORBIT)) cerr << "ORBIT Uninitialized" << endl;
  if(!(status & POSITION)) cerr << "POSITION Uninitialized" << endl;
  if(!(status & SOURCEDIR)) cerr << "SOURCEDIR Uninitialized" << endl;
  if(!(status & PARALLAX)) cerr << "PARALLAX Uninitialized" << endl;
  if(!(status & PARAMS)) cerr << "PARAMS Uninitialized" << endl;
  if(!(status & EPOCHS)) cerr << "EPOCHS Uninitialized" << endl;
  if(!(status & PLLXINITIALIZED)) cerr << "PLLXINITIALIZED Uninitialized" << endl;
}

void parallax::compute_NEshifts()
{
  if((status&(allready-PLLXINITIALIZED)) != (allready-PLLXINITIALIZED))
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Some element of the parallax computation is uninitialized:" << endl;
      print_uninit();
      exit(1);
    }

  vector<double> x, xp;
  NEshift = vector<vector<double> >(epochs.size(),vector<double>(2,0));

  for(int i=0;i<int(epochs.size());i++)
    {
      //compute the shift relative to the reference frame
      x = vector<double>(3,0);
      //orbit->viewfrom(epochs[i],a,d,&x);
      for(int j=0;j<int(orbit->size());j++)
	{
	  (*orbit)[j].viewfrom(epochs[i],a,d,&xp);
	  for(int k=0;k<int(x.size());k++) x[k] += xp[k];
	  //figure out how to print out the x,y,z coordinates in the ecliptic plane here
	}
      NEshift[i][0] = x[0] - xref[0] - (epochs[i]-tref)*vref[0];
      NEshift[i][1] = x[1] - xref[1] - (epochs[i]-tref)*vref[1];
    }
}

void parallax::compute_tushifts()
{
  if(status!=allready)
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Some element of the parallax computation is uninitialized" << endl;
      exit(1);
    }

  tushift = vector<vector<double> >(epochs.size(),vector<double>(2,0));

  double cs= cos(phi_pi); 
  double sn= sin(phi_pi);

  for(int i=0;i<int(epochs.size());i++)
    {
      //Convert the shift in the observer plane to a shift in the 
      //source position

      tushift[i][0] = -piE * ( NEshift[i][0]*cs + NEshift[i][1]*sn);
      tushift[i][1] = -piE * (-NEshift[i][0]*sn + NEshift[i][1]*cs);
    }
}


void parallax::setup_reference_frame() //private
{
  //Assumes that all setup needed has been done with the orbits
  int exitnow=0;

  if(!(status&REFFRAME))
    {
      cerr << __FILE__ << ": Error: Reference frame not set" << endl;
      exitnow=1;
    }

  if(!(status&ORBIT))
    {
      cerr << __FILE__ << ": Error: Orbit not set" << endl;
      exitnow=1;
    }

  if(!(status&POSITION))
    {
      cerr << __FILE__ << ": Error: Event position not set" << endl;
      exitnow=1;
    }

  if(exitnow) exit(1);

  //Set up the reference frame
  //oref->viewfrom(tref,a,d,&xref);
  //oref->velocity(tref,a,d,&vref);
  //oref->acceleration(tref,a,d,&aref);
  
  //reset the orbits
  xref=vref=aref=vector<double>(3,0);

  //add on purturbations
  vector<double> xp, vp, ap;
  for(int i=0;i<int(oref->size());i++)
    {
      (*oref)[i].viewfrom(tref,a,d,&xp);
      (*oref)[i].velocity(tref,a,d,&vp);
      (*oref)[i].acceleration(tref,a,d,&ap);
      for(int j=0;j<int(xp.size());j++) xref[j]+=xp[j];
      for(int j=0;j<int(vp.size());j++) vref[j]+=vp[j];
      for(int j=0;j<int(ap.size());j++) aref[j]+=ap[j];
    }
}

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
  //(velocity) in target frame
  
  //wr = (piEN,piEE)/tE is the lens velocity in Einstein radii per day
  //wref = apparent motion of the reference frame in the same units, as viewed 
  //from the sun

  wref[0] = vref[0]*piE;
  wref[1] = vref[1]*piE;

  //Calculate proper motion components in NE and llrp
  phi_llN = atan2(-aref[1],-aref[0]);

  //set_pm_lb(mul_,mub_)      progress[LSOURCEDIR]=1;
  //set_pm_ad(mua_,mud_)      progress[LSOURCEDIR]=1;
  //set_piEpp(piEll_,piErp_)  progress[LSOURCEDIR]=2;
  //set_piENE(piEN_,piEE_)    progress[LSOURCEDIR]=3;

  //cerr << progress[LPARAMS] << " " << progress[LSOURCEDIR] << endl;

  if(progress[LPARAMS]==1 && progress[LSOURCEDIR]==1)
    {
      wh[0] = mud/tE_h;
      wh[1] = mua/tE_h;
      wr[0] = wh[0] - wref[0];
      wr[1] = wh[1] - wref[1];
      tE_r = 1.0/qAdd(wr[0],wr[1]);
      piEN = wr[0]*tE_r;
      piEE = wr[1]*tE_r;
      piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
      piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);
      progress[LSOURCEDIR]=2;
    }
  else if(progress[LPARAMS]==2 && progress[LSOURCEDIR]>1)
    {
      if(progress[LSOURCEDIR]==2)
	{
	  piEll = piEN*cos(-phi_llN) - piEE*sin(-phi_llN);
	  piErp = piEN*sin(-phi_llN) + piEE*cos(-phi_llN);
	}
      else if(progress[LSOURCEDIR]==3)
	{
	  piEN = piEll*cos(phi_llN) - piErp*sin(phi_llN);
	  piEE = piEll*sin(phi_llN) + piErp*cos(phi_llN);
	}
      wr[0] = piEN/tE_r;
      wr[1] = piEE/tE_r;
      wh[0] = wr[0] + wref[0];
      wh[1] = wr[1] + wref[1];
      tE_h = 1.0/qAdd(wh[0],wh[1]);
      mud = wh[0]*tE_h;
      mua = wh[1]*tE_h;
      c.muad2lb(a,d,mua,mud,&mul,&mub);
    }
  else
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Must specify parameters and source directions as heliocentric or reference-frame-centric" << endl;
    }

  phi_pi = atan2(piEE,piEN);
  
}

int parallax::set_reference(double tref_, vector<orbitalElements>* oref_)
{
  tref=tref_;
  oref = oref_;

  status |= REFFRAME;
  return status;
}

int parallax::set_orbit(vector<orbitalElements>* orbit_)
{
  orbit = orbit_;
  status |= ORBIT;
  return status;
}

int parallax::set_radec(double ra, double dec, int deghr)
{
  //deghr = 0 - both radians
  //      = 1 - both degrees
  //      = 2 - hr and degrees

  if(deghr>2)
    {
      cerr << __FILE__ << ": "<< __FUNCTION__ << ": Error setting position - bad angle unit code (" << deghr << ")" << endl;
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

  status |= POSITION;
  return status;
}

int parallax::set_radec(vector<double> ra, vector<double> dec)
{
  //deghr = 0 - both radians
  //      = 1 - both degrees
  //      = 2 - hr and degrees

  if(int(ra.size())!=3||int(dec.size())!=3)
    {
      cerr << __FILE__ << ": "<< __FUNCTION__ << ": Error setting position - input sexagesimal vectors (ra_hr,ra_min,ra_sec), (dec_deg,dec_min,dec_sec)" << endl;
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

  status |= POSITION;
  return status;
}

int parallax::set_lb(double l_, double b_, int deg)
{
  //deg = 0 - both degrees
  //    = 1 - both radians

  if(deg>1)
    {
      cerr << __FILE__ << ": "<< __FUNCTION__ << ": Error setting position - bad angle unit code (" << deg << ")" << endl;
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
  
  status |= POSITION;
  return status;
}

int parallax::set_pm_lb(double mul_, double mub_)
{
  //set the heliocentric proper motion of the lens relative to the source in 
  //Galactic coordinates
  //assume mu_l = mu_l cos b ALWAYS
  
  if(!(status & POSITION))
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << "Must set position in order to do proper motion unit conversions" << endl;
      exit(1);
    }
    
  double mu = qAdd(mul_,mub_);
  mul = mul_/mu;
  mub = mub_/mu;
  c.mulb2ad(l,b,mul,mub,&mua,&mud);

  progress[LSOURCEDIR]=1;
  status |= SOURCEDIR;
  return status;
  
}

int parallax::set_pm_ad(double mua_, double mud_)
{
  //set the heliocentric proper motion of the lens relative to the source in 
  //celestoal coordinates
  //assume mu_a = mu_a cos d ALWAYS
  
  if(!(status & POSITION))
    {
      cerr << __FILE__ << ": " << __FUNCTION__ << "Must set position in order to do proper motion unit conversions" << endl;
      exit(1);
    }
    
  double mu = qAdd(mua_,mud_);
  mua = mua_/mu;
  mud = mud_/mu;
  c.muad2lb(a,d,mua,mud,&mul,&mub);

  progress[LSOURCEDIR]=1;
  status |= SOURCEDIR;
  return status;
  
}

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

int parallax::set_piE(double piE_)
{
  piE = piE_;

  status |= PARALLAX;
  return status;
}

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

int parallax::set_tE_h(double tE)
{
  tE_h = tE;

  progress[LPARAMS]=1;
  status |= PARAMS;
  return status;
}

int parallax::set_tE_r(double tE)
{
  tE_r = tE;

  progress[LPARAMS]=2;
  status |= PARAMS;
  return status;
}

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
