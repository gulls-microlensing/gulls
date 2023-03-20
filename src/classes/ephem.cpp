#include <cmath>
#include <vector>

#include "ephem.h"
#include "constants.h"
#include "integerPowers.h"

using namespace std;

//set up some common orbits

void orbitalElements::blank()
{
  //a stationary object
    
  epoch=J2000;
  a0=e0=I0=L0=w0=O0=0;
  da=de=dI=dL=dw=dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::earthmoonbary()
{
  epoch=J2000;
  a0=1.00000261; da=0.0000000562;
  e0=0.01671123; de=-0.0000004392;
  I0=-0.00001531; dI=-0.0001294668;
  L0=100.46457166; dL=359.9937244981;
  w0=102.93768193; dw=0.0032327364;
  O0=0.0; dO=0.0;
  b=0; c=0; s=0; f=0;
      
  storerad();
};

void orbitalElements::earth()
{
  //Orbit of the earth center about the earth moon barycenter
  //Taken from JPL horizons for 01 Jan 2022
  epoch=2459580.5;
  a0=3.129259146938130E-05; da=0;
  e0=7.097579690698613E-02; de=0;
  I0=5.182392994073342E+00; dI=0;
  L0=3.458150871901001E+02; dL=1.314938028470777E+01*365.25;
  w0=3.064812761639257E+01; dw=0;
  O0=6.093679521123980E+01; dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::earthl2()
{
  //sun/earth L2 - identical to Earth, except for a slightly larger a (changing a does not change the period)
  earthmoonbary();
  double shift = 1 + pow(Mearth/3.0,1.0/3.0);
  a0*=shift;
  da*=shift;
}

void orbitalElements::geosynch(double inc, double phase)
{
  //roughly corresponds to the GOES orbit
  //phase is between 0 and 1
  epoch=2459581.0;
  a0=42168.21876/aukm;    da=0;
  e0=0.0002;              de=0;
  I0=23.43927944+inc;     dI=0;
  L0=phase*360.0;         dL=131850.0538;
  w0=120.0;               dw=0;
  O0=30;                  dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::lissajous(double phase)
{
  //bad approximation of a Lissajous orbit about L2
  //Based on orbit of WMAP taken from SNAP-TECH-04010
  //phase is between 0 and 1
  epoch=2459581.0;
  a0=276000/aukm;         da=0;
  e0=0;                   de=0;
  I0=73;                  dI=0;
  L0=phase*360.0;         dL=738.7078652; //178 day period
  w0=0;                   dw=0;
  O0=90;                  dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::lissajousxy(double phase)
{
  //rough approximation of a Lissajous orbit about L2 in the xy plane, combine
  //with lissajousz
  //Based on orbit of WMAP taken from SNAP-TECH-04010
  //phase is between 0 and 1
  epoch=2459581.0;
  a0=82000/aukm;          da=0;
  e0=0;                   de=0;
  I0=0;                   dI=0;
  L0=phase*360.0;         dL=740.3716216; //177.6 day period
  w0=0;                   dw=0;
  O0=0;                   dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::lissajousz(double phase)
{
  //rough approximation of a Lissajous orbit about L2 in the z plane, combine
  //with lissajousxy
  //Based on orbit of WMAP taken from SNAP-TECH-04010
  //phase is between 0 and 1
  epoch=2459581.0;
  a0=264000/aukm;         da=0;
  e0=0;                   de=0;
  I0=90;                  dI=0;
  L0=phase*360.0;         dL=740.3716216; //177.6 day period
  w0=0;                   dw=0;
  O0=90;                  dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}

void orbitalElements::mars()
{
  epoch=J2000;
  a0=1.52371034;   da=0.0000001847;
  e0=0.09339410;   de=0.0000007882; 
  I0=1.84969142;   dI=-0.0000813131;
  L0=-4.55343205;  dL=191.4030268499;
  w0=-23.94362959; dw=0.0044441088;
  O0=49.55953891;  dO=-0.0029257343;
  b=0; c=0; s=0; f=0;

  storerad();
};

void orbitalElements::jupiter()
{
  epoch=J2000;
  a0=5.20288700;   da=-0.0000011607;
  e0=0.04838624;   de=-0.0000013253; 
  I0=1.30439695;   dI=-0.0000183714; 
  L0=34.39644051;  dL=30.3474612775;
  w0=14.72847983;  dw=0.0021252668;
  O0=100.47390909; dO=0.0020469106;
  b=0; c=0; s=0; f=0;

  storerad();
};

void orbitalElements::jwst(double phase)
{
  //bad approximation of JWST's orbit about L2
  //Based on orbit of WMAP taken from SNAP-TECH-04010
  //phase is between 0 and 1
  epoch=2459581.0;
  a0=0.8e6/aukm;          da=0;
  e0=0;                   de=0;
  I0=73;                  dI=0;
  L0=phase*360.0;         dL=738.7078652; //178 day period
  w0=0;                   dw=0;
  O0=90;                  dO=0;
  b=0; c=0; s=0; f=0;

  storerad();
}


//Now to the calculations

//set the julian date, and work out the values of the elements
void orbitalElements::settime(double jd)
{
  t = jd;
  T = (t-epoch)/365.25;

  a = a0 + T*da;
  e = e0 + T*de;
  I = I0 + T*dI;
  L = L0 + T*dL;
  w = w0 + T*dw;
  O = O0 + T*dO;
};


//compute the mean anomaly and arg of perihelion
void orbitalElements::meananomaly()
{
  W = w - O;
  //mean anomaly - dL handles its change with time
  M = L - w + b*T*T + c*cos(f*T) + s*sin(f*T);

  //put into the range -180<=M<=180
  //first get into 0,360
  if(M>twoPi) M = fmod(M,twoPi);
  else if(M<-pi) M = fmod(M-floor(M/twoPi)*twoPi,twoPi);

  //now get to -180,180
  if(M>pi) M-=twoPi;
}
  
//compute the eccentric anomaly
void orbitalElements::eccentricanomaly()
{
  //iterative solution of Kepler's equations
  double dM, dE;

  E = M + e * sin(M);

  do
    {
      dM = M - (E - e*sin(E));
      dE = dM/(1 - e*cos(E));
      E += dE;
    } while(abs(dE)> tol);
}    

//compute the coordinates in the plane of the orbit
void orbitalElements::heliocoords()
{
  xh = a*(cos(E)-e);
  yh = a*sqrt(1-e*e)*sin(E);
  zh = 0;
}

//compute coordinates in the J2000 ecliptic plane
void orbitalElements::eclcoords()
{
  double cW = cos(W); double sW = sin(W);
  double cO = cos(O); double sO = sin(O);
  double cI = cos(I); double sI = sin(I);
  
  xecl = (cW*cO - sW*sO*cI) * xh + (-sW*cO - cW*sO*cI) * yh;
  yecl = (cW*sO + sW*cO*cI) * xh + (-sW*sO + cW*cO*cI) * yh;
  zecl = (sW*sI) * xh + (cW*sI) * yh;
}

//compute coordinates in the equatorial ICRS
void orbitalElements::eqcoords()
{
  xeq = xecl;
  yeq = cosob*yecl - sinob*zecl;
  zeq = sinob*yecl + cosob*zecl;
}

void orbitalElements::compute()
{
  meananomaly();
  eccentricanomaly();
  heliocoords();
  eclcoords();
  eqcoords();
}

void orbitalElements::perturb(double t, vector<orbitalElements*> perturbers)
{
  settime(t);
  compute();

  for(int i=0;i<int(perturbers.size());i++)
    {
      perturbers[i]->settime(t);
      perturbers[i]->compute();
      perturb(perturbers[i]);
    }
  
}

void orbitalElements::perturb(orbitalElements* perturbation)
{
  //perturb the 3D position of the orbit by the 3D position stored in the passed elements
  //i.e. run compute in both elements, then call perturb

  //heliocentric coordinates in the plane of orbit
  xh += perturbation->xh;
  yh += perturbation->yh;
  zh += perturbation->zh;

  //coordinates in the J2000 ecliptic plane
  xecl += perturbation->xecl;
  yecl += perturbation->yecl;
  zecl += perturbation->zecl;

  //coordinates in J2000 equatorial plane (ICRF)
  xeq += perturbation->xeq;
  yeq += perturbation->yeq;
  zeq += perturbation->zeq;

}

void orbitalElements::viewfrom(double t, double ra, double dec, vector<double>* view)
{
  //one stop shop for computing position at time t, then viewing from the sun 
  //in the ra,dec (radians) direction

  settime(t);
  compute();
  viewfrom(ra, dec, view);
}

void orbitalElements::viewfrom(double ra, double dec, vector<double>* view)
{
  //view the orbit from the sun, looking in the direction of ra,dec, with 
  //the sun at the origin - units in AU, ra and dec in radians, normal (n) axis 
  //points from given ra, dec to sun

  //Start off working in the cartesian celestial frame where x points to the vernal equinox, z to the NCP, y in the direction of increasing RA

  double cr = cos(ra); double cd = cos(dec);
  double sr = sin(ra); double sd = sin(dec);

  //first work out the ra,dec unit vector (normal) nhat
  //This is a cartesian unit vector pointing from the sun to the ra,dec point (i.e. the oposite direction to the viewing direction
  vector<double> nhat(3,0);
  vector<double> Nhat(3,0);
  vector<double> Ehat(3,0);
  vector<vector<double> > T(3, vector<double>(3,0));
  nhat[0]=cr*cd;          
  nhat[1]=sr*cd;
  nhat[2]=sd; 

  //In the sky frame:
  //the east vector on the sky will be perpendicular to both the normal and NCP vectors
  //i.e. East = zhat x nhat - will need to normalize as these two are not perpendicular
  double length=qAdd(nhat[0],nhat[1]);
  if(length>0)
    {
      Ehat[0]=-nhat[1]/length;
      Ehat[1]=nhat[0]/length;
    }
  else
    {
      //arbitrary definition in the undefined case or norm || north
      //maybe it shouldn't be arbitrary, but this shouldn't come up
      Ehat[0]=1;
      Ehat[1]=0;
    }
  Ehat[2]=0;

  //the North vector on the sky will be perpendicular to both east and the normal
  //i.e. North = norm x East, both are unit and perpendicular, so no need to normalize  
  Nhat[0]=-nhat[2]*Ehat[1];
  Nhat[1]=nhat[2]*Ehat[0];
  Nhat[2]=nhat[0]*Ehat[1]-nhat[1]*Ehat[0];

  //The N,E,n coordinate system is now a rhs with the barycenter as the origin. E and N are the cardinal directions on the sky, and n is the distance. Positive n implies the object is closer to the star at ra,dec than the barycenter is.

  //Find the position of the object in the ENn frame
  vector<double> obj_NEn(3);
  obj_NEn[0] = Nhat[0]*xeq + Nhat[1]*yeq + Nhat[2]*zeq;
  obj_NEn[1] = Ehat[0]*xeq + Ehat[1]*yeq + Ehat[2]*zeq;
  obj_NEn[2] = nhat[0]*xeq + nhat[1]*yeq + nhat[2]*zeq;

  view->resize(3);
  (*view) = obj_NEn;
}


void orbitalElements::viewfrom_old(double ra, double dec, vector<double>* view)
{
  //view the orbit from a given ra,dec, with the sun at the origin - units in AU
  //ra and dec in radians, z axis points from given ra, dec to sun

  double cr = cos(ra); double cd = cos(dec);
  double sr = sin(ra); double sd = sin(dec);

  //first work out the ra,dec unit vector (normal)
  //This is a cartesian unit vector pointing from the sun to the ra,dec point
  vector<double> norm(3);
  vector<double> north(3);
  vector<double> east(3);
  norm[0]=cr*cd;          
  norm[1]=sr*cd;
  norm[2]=sd;  //currently in the cartesian frame fixed to celestial sphere

  //In the cartesian frame:
  //the North vector will always point to 0,0,1
  //the East vector will always point to 1,0,0
  //In the sky frame:
  //the east vector on the sky will be perpendicular to both the normal and north vectors
  //i.e. east = North x norm - will need to normalize as these two are not perpendicular
  double length=qAdd(norm[0],norm[1]);
  if(length>0)
    {
      east[0]=-norm[1]/length;
      east[1]=norm[0]/length;
    }
  else
    {
      //arbitrary definition in the undefined case or norm || north
      //maybe it shouldn't be arbitrary, but this shouldn't come up
      east[0]=1;
      east[1]=0;
    }
  east[2]=0;

  //the north vector on the sky will be perpendicular to both east and the normal
  //i.e. north = norm x east, both are unit and perpendicular, so no need to normalize  
  north[0]=-norm[2]*east[1];
  north[1]=norm[2]*east[0];
  north[2]=norm[0]*east[1]-norm[1]*east[0];

  //Now, project the position of the planet onto the sky frame
  view->resize(2);
  (*view)[0] = xeq*east[0] + yeq*east[1] + zeq*east[2];
  (*view)[1] = xeq*north[0] + yeq*north[1] + zeq*north[2];
}

void orbitalElements::velocity(double t, double ra, double dec, vector<double>* v, double tstep)
{
  vector<double> x1(3),x2(3);
  v->resize(3);

  //tstep is fraction of a period
  tstep *= 365.25*360.0/dL; //in days  
  
  viewfrom(t-tstep,ra,dec,&x1);
  viewfrom(t+tstep,ra,dec,&x2);

  (*v)[0] = ((x2[0]-x1[0])/(2*tstep)); //*AUday2kms;
  (*v)[1] = ((x2[1]-x1[1])/(2*tstep)); //*AUday2kms;
  (*v)[2] = ((x2[2]-x1[2])/(2*tstep))*AUday2kms;
}

void orbitalElements::acceleration(double t, double ra, double dec, vector<double>* acc, double tstep)
{
  vector<double> x1(3),x2(3),x3(3);
  acc->resize(3);

  //tstep is fraction of a period
  tstep *= 365.25*360.0/dL; //in days
  double tstep2 = tstep*tstep;
  
  viewfrom(t-tstep,ra,dec,&x1);
  viewfrom(t,ra,dec,&x2);
  viewfrom(t+tstep,ra,dec,&x3);

  (*acc)[0] = ((x3[0]+x1[0]-2*x2[0])/(tstep2)); //*AUday2kms; in AU/d^2
  (*acc)[1] = ((x3[1]+x1[1]-2*x2[1])/(tstep2)); //*AUday2kms;
  (*acc)[2] = ((x3[2]+x1[2]-2*x2[2])/(tstep2))*AUday2kms/day2s; //in km s^-2
}

void orbitalElements::eqpropermotion(double t, double ra, double dec, double dist, vector<double>* pm, double tstep)
{
  //dist in kpc
  //Proper motion of the planet as viewed from the sun towards the direction specified. Returns pm_ra*cos(dec),pm_dec,RV is mas yr-1,km s-1

  vector<double> v(2);

  double scale=365.25/AUday2kms/dist;

  velocity(t, ra, dec, &v, tstep);

  pm->resize(2);
  (*pm)[0] = scale*v[0]*cos(dec);
  (*pm)[1] = scale*v[1];
  (*pm)[2] = v[2];
  
}

void orbitalElements::galpropermotion(double t, double ll, double bb, double dist, vector<double>* pmgal, double tstep)
{
  //dist in kpc
  //Proper motion of the planet as viewed from the sun towards the direction specified. Returns pm_ra*cos(dec),pm_dec is mas yr-1

  vector<double> gal(2);
  vector<double> eq(2);
  vector<double> pmeq(2);

  gal[0]=ll; gal[1]=bb;  
  gal2eq(gal, &eq);
  
  eqpropermotion(t, eq[0], eq[1], dist, &pmeq, tstep);
  pmeq2gal(eq,pmeq,pmgal);
}

void orbitalElements::pmeq2gal(vector<double> eq, vector<double> pmeq, vector<double>* pmgal)
{
  double sinp, cosp;
  double cosb, sinb;
  double sind, cosd;

  vector<double> gal(2);

  eq2gal(eq, &gal);
  sind = sin(eq[1]);
  cosd = cos(eq[1]);
  cosb = cos(gal[1]);
  sinb = sin(gal[1]);

  cosp = (sindngp - sind*sinb)/(cosd*cosb);
  sinp = sin(eq[0] - rangp);

  pmgal->resize(2);
  (*pmgal)[0] =  cosp*pmeq[0] + sinp*pmeq[1];
  (*pmgal)[1] = -sinp*pmeq[0] + cosp*pmeq[1];

}

void orbitalElements::pmgal2eq(vector<double> gal, vector<double> pmgal, vector<double>* pmeq)
{
  double sinp, cosp;
  double cosb, sinb;
  double sind, cosd;

  vector<double> eq(2);

  gal2eq(gal, &eq);
  cosb = cos(gal[1]);
  sinb = sin(gal[1]);
  sind = sin(eq[1]);
  cosd = cos(eq[1]);

  cosp = (sindngp-sind*sinb)/(cosd*cosb);
  sinp = sin(eq[0]-rangp)*cosdngp/cosb;

  pmeq->resize(2);
  (*pmeq)[0] = cosp*pmgal[0] + sinp*pmgal[1];
  (*pmeq)[1] = -sinp*pmgal[0] + cosp*pmgal[1];
}

void orbitalElements::eq2gal(vector<double> eq, vector<double>* gal)
{
  double sinb, cosb, sinlml0, coslml0;
  double sind, cosd, sinama0, cosama0;

  double ra=eq[0];
  double dec=eq[1];

  sind=sin(dec); cosd=cos(dec);
  sinama0=sin(ra-ra0); cosama0=cos(ra-ra0);

  sinb = sind*sindngp - cosd*cosdngp*sinama0;
  cosb = sqrt(1-sinb*sinb);
  coslml0 = cosama0*cosd/cosb;
  sinlml0 = (sind*cosdngp + cosd*sindngp*sinama0)/cosb;
  
  gal->resize(2);
  (*gal)[0] = ll0 + atan2(sinlml0,coslml0);
  (*gal)[1] = asin(sinb);
}
 
void orbitalElements::gal2eq(vector<double> gal, vector<double>* eq)
{
  double sinb, cosb, sinlml0, coslml0;
  double sind, cosd, sinama0, cosama0;

  double ll=gal[0];
  double bb=gal[1];

  sinb=sin(bb); cosb=cos(bb);
  sinlml0=sin(ll-ll0); coslml0=cos(ll-ll0);

  sind = sinb*sindngp + cosb*cosdngp*sinlml0;
  cosd = sqrt(1.0-sind*sind);
  cosama0 = coslml0*cosb/cosd;
  sinama0 = (-sinb*cosdngp + cosb*sindngp*sinlml0)/cosd;  

  eq->resize(2);
  (*eq)[0] = ra0 + atan2(sinama0,cosama0);
  (*eq)[1] = asin(sind);
  
  
}
