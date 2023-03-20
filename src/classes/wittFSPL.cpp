#include<cmath>
#include<complex>

#include<gsl/gsl_sf_ellint.h>
#include<gsl/gsl_mode.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>


#include"constants.h"
#include"integerPowers.h"
#include"wittFSPL.h"
#include"singleLens.h"
#include"src_cld.h"

////////////////////////////////////////////////////////////
//
//    Finite source point lens magnification
//      taken from Witt & Mao, 1994, ApJ, 430, 505
//      http://adsabs.harvard.edu/abs/1994ApJ...430..505W
//
///////////////////////////////////////////////////////////


double wittFSMagnification(double zs0, double rs)
{
  //cout << "entering witt "; cout.flush();
  double rs2=sqr(rs);
  double sub=zs0-rs;

  if(abs(sub)<1.0e-7)   //when zs0=rs
    {
      return (1.0/pi) * ( 2.0/rs
			  + (1.0+rs2)*(piO2 + asin((rs2-1.0)/(rs2+1.0)))/rs2);

    }
  else   //the standard case
    {
      double add=zs0+rs;
      double q=qAdd(2.0,sub);
  

      double n=4.0*rs*zs0/sqr(add); 
      double k=sqrt(4.0*n)/q; 

      double np=-n;   //convert to the correct convention for GSL

      gsl_mode_t mode=GSL_PREC_DOUBLE;

      //E(pi/2,k) = Ecomp(k)
      double E = gsl_sf_ellint_Ecomp(k,mode);
      
      //F(pi/2,k) = Kcomp(k)
      double F = gsl_sf_ellint_Kcomp(k,mode);
      
      double P = gsl_sf_ellint_Pcomp(k,np,mode);

      double mm=(1.0/(twoPi*rs2)) * ( 0.5*E*add*q
				      - F*sub*(4.0+0.5*(sqr(zs0)-rs2))/q
				      + P*2.0*sqr(sub)*(1.0+rs2)/(q*add) 
				      );

      //double mu1 = 0.5 + mm;
      //double mu2 = -0.5 + mm;      
      //return abs(mu1)+abs(mu2);
      if(mm<-0.5) return -2.0*mm;
      else if(mm>=0.5) return 2.0*mm;
      else return 1.0;
   }
}



void dmududr_witt(double u, double rs, double* dmudu, double* dmudr, double delta)
{
  //for |u|>10rs, uses the Hexadecapole approximation, for u<10 rs, does a numerical derivative of the Witt & Mao magnification

  static const double ld2=1.0;
  static const double ld4=1.0;

  double logu = log10(abs(u));
  double logrs = log10(rs);
  
  double up = u*(1+delta);
  double um = u*(1-delta);
  *dmudu = 0.5*(
		wittFSMagnification(abs(up),rs)
		-wittFSMagnification(abs(um),rs)
		)/(abs(u)*delta);
  *dmudr = 0.5*(wittFSMagnification(abs(u),rs*(1+delta))
		-wittFSMagnification(abs(u),rs*(1-delta))
		)/(rs*delta); 

  //Stop numerical errors
  if(logu-logrs>0.5 && logu-logrs/3.0>1.5)
    {
      //Use the PSPL derivative for dmudu
      *dmudu = pspl_dmudu(u);
    }

  if(logu-logrs>0.5 && logu-2.0*logrs/3.0>1)
    {
      //Use the PSPL derivative for dmudu
      static const double cc = -1.94299; //         +/- 0.002516     (0.1295%)
      //static const double kk = 2.82854; //          +/- 0.004476     (0.1582%)
      static const double logkk = 0.451562325;

      *dmudr = pow(10,logrs + cc + (logu<logkk?-3:-5.85)*(logu-logkk));
    }
 

  /*  if(abs(u)>10*rs && rs<1e-2)
    {
      vector<double> dAiidu(13);
      vector<double> dAiidr(13);
      vector<double> Aii(13);
      for(int i=0;i<13;i++)
	{
	  double uj = qAdd(u+gphOffX[i]*rs, gphOffY[i]*rs);
	  double dAduj = pspl_dmudu(uj);
	  double dujdu = (u+gphOffX[i]*rs)/uj;
	  double dujdr = ((u+gphOffX[i]*rs)*gphOffX[i]+sqr(gphOffY[i])*rs)/uj;
	  Aii[i] = pacAmp(uj);
	  if(i>0)
	    {
	      dAiidu[i] = dAduj*dujdu;
	      dAiidr[i] = dAduj*dujdr;
	    }
	  else
	    {
	      dAiidu[0] = dAduj;
	      dAiidr[0] = 0;
	    }
	}
      double dAhdu = 0.25*(dAiidu[1]+dAiidu[2]+dAiidu[3]+dAiidu[4])
	-dAiidu[0];
      double dApdu = 0.25*(dAiidu[5]+dAiidu[6]+dAiidu[7]+dAiidu[8])
	-dAiidu[0];
      double dAxdu = 0.25*(dAiidu[9]+dAiidu[10]+dAiidu[11]+dAiidu[12])
	-dAiidu[0];
      
      double dAhdr = 0.25*(dAiidr[1]+dAiidr[2]+dAiidr[3]+dAiidr[4])
	-dAiidr[0];
      
      double dApdr = 0.25*(dAiidr[5]+dAiidr[6]+dAiidr[7]+dAiidr[8])
	-dAiidr[0];
      double dAxdr = 0.25*(dAiidr[9]+dAiidr[10]+dAiidr[11]+dAiidr[12])
	-dAiidr[0];
      
      double dA2rdu = (16*dAhdu-dApdu)/3.0;
      double dA2rdr = (16*dAhdr-dApdr)/3.0;

      double dA4rdu = 0.5*(dApdu+dAxdu) - dA2rdu;
      double dA4rdr = 0.5*(dApdr+dAxdr) - dA2rdr;

      *dmudu = dAiidu[0] + 0.5*ld2*dA2rdu + ld4/3.0*dA4rdu;
      *dmudr = dAiidr[0] + 0.5*ld2*dA2rdr + ld4/3.0*dA4rdr;
    }
  else
    {
      if(u>20&&u/rs>20)
	{
	  *dmudu = pspl_dmudu(u);
	  *dmudr = 0;
	  //double up = u*(1+delta);
	  //double um = u*(1-delta);
	  // *dmudu = 0.5*(wittFSMagnification(abs(up),rs)-wittFSMagnification(abs(um),rs))/(abs(u)*delta);
	  // *dmudr = 0.5*(wittFSMagnification(abs(u),rs*(1+delta))-wittFSMagnification(abs(u),rs*(1-delta)))/(rs*delta); 
	}
      else
	{
	  double up = u*(1+delta);
	  double um = u*(1-delta);
	  *dmudu = 0.5*(
			wittFSMagnification(abs(up),rs)
			-wittFSMagnification(abs(um),rs)
			)/(abs(u)*delta);
	  *dmudr = 0.5*(wittFSMagnification(abs(u),rs*(1+delta))
			-wittFSMagnification(abs(u),rs*(1-delta))
			)/(rs*delta); 
	}
    }
  */
}


////////////////////////////////////////////////////////////
//
//    Finite source point lens magnification
//      taken from Lee et al., 2010, ApJ, 695, 200
//      http://adsabs.harvard.edu/abs/2010ApJ...695..200L
//
///////////////////////////////////////////////////////////


double leeFSMagnification(double u, double r)
{

  double lee_f(double phi, double u, double r);
  double lee_fg(double phi, double u, double r, double asinpu);
  //cout << "entering witt "; cout.flush();
  double r2=sqr(r);
  double sub=u-r;
  double add=u+r;
  double qp=qAdd(add,2);
  double qm=qAdd(sub,2);

  double mu;
  double mu0 = (add*qp-sub*qm);
  double mu1=0, mu2=0;
  int N=100; //N must be even
  double pioN=pi/N;

  if(u<=r)
    {
      for(int k=1;k<=N-1;k++)
	{
	  mu1+=lee_f(k*pioN,u,r);
	  mu2+=lee_f((k-0.5)*pioN,u,r);
	}
      mu2+=lee_f((N-0.5)*pioN,u,r);
      mu = (mu0+2*mu1+4*mu2)/(6.0*r2*N);
    }
  else
    {
      double asinpu=asin(r/u);
      double n=N/2;

      for(int k=1;k<=n-1;k++)
	{
	  mu1+=lee_fg(2*k*asinpu/N,u,r,asinpu);
	  mu2+=lee_fg((2*k-1)*asinpu/N,u,r,asinpu);
	}
      mu2+=lee_fg((2*n-1)*asinpu/N,u,r,asinpu);
      mu = asinpu*(mu0+2*mu1+4*mu2)/(3.0*pi*r2*N);
    }
  
  return mu;
}

//f when u<=r
double lee_f(double theta, double u, double r)
{
  double u2, uct;

  uct = u*cos(theta);
  u2 = uct + sqrt(r*r-u*u+uct*uct);

  return u2*qAdd(u2,2);  
}

//f when u>r
double lee_fg(double theta, double u, double r, double asinpu)
{
  double u1, u2, uct;

  if(theta > asinpu) return 0;

  uct = u*cos(theta);
  u1 = uct - sqrt(r*r-u*u+uct*uct);
  u2 = uct + sqrt(r*r-u*u+uct*uct);
  return u2*qAdd(u2,2) - u1*qAdd(u1,2);  
}



double fsplMagnification(double zs, src_cld* src) 
{
  int i;
  struct wfspl_parameters p;
  double drs, rs0;
  double area=0.0;
  double limb;
  double u1, u2;
  double R;
  double mu;
  const int nr=100;

  /* GSL variables */
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  
  double epsAbsolute=1.0e-5;
  double epsRelative=0.0;//1.0e-5;
  double result, error;

  gsl_function F;
  int status;

  double Agph, A2;

  p.zs = zs;
  p.Kx = 0.0;
  p.Ky = 0.0;

  drs = src->rs/nr;
  rs0 = -drs * 0.5;


  u1 = src->ldGamma; u2=0.0; /* linear LD  profile */ 


   /*ATTEMPT HERE TO SPEED UP FITTER BY SETTING MAGNIFICATION 
      FOR STARS > 10*rho AWAY FROM LOS TO LENS TO PS solution AND NOT CALLING
      INTEGRATION */

  Agph = gphmag(zs,src,A2);
  static int swap=0;

  /*  if(zs>10.0*rs)
    {
      mu = (zs*zs + 2)/(zs*sqrt(zs*zs + 4));
      gsl_integration_workspace_free (w);
      }*/
  if(abs(A2)<epsAbsolute)
    {
      gsl_integration_workspace_free (w);
      mu=Agph;
      if(swap!=1) cout << zs << endl;
      swap=1;
    }
  else
    {
      if(swap!=2) cout << zs << endl;
      swap=2;

      mu = area = 0.0;

      for (i=0; i<nr; i++) 
	{
	  rs0 += drs;
	  p.rs = rs0;
	  
	  F.params = &p;
	  
	  R = rs0/src->rs;
	  limb = 1.0-u1-u2+u1*sqrt(1.-R*R)+u2*(1.0-R*R);  
	  
   
	  gsl_set_error_handler_off();
	  /* integrate to obtain magnification */
	  F.function = &muFunc;
 
	  status = gsl_integration_qags (&F, 0.0, twoPi, epsAbsolute, epsRelative, 10000, w, &result, &error);    if(status) {cerr << "wittFSPL.cpp: fsplMagnification: Error caught: " << status << endl; exit(1);}
    

	  mu += result * limb;

 
	  area += twoPi*rs0 * limb;
	}

      gsl_integration_workspace_free (w);

 
      mu /= area;
      /* printf("%f %f %f\n",zs,*mu,(zs*zs + 2)/(zs*sqrt(zs*zs + 4))); */
    }

  return mu;
} 



/* function to be integrated for magnification */
double muFunc(double phi, void *parameter)
{
  struct wfspl_parameters *p = (struct wfspl_parameters *) parameter;
  double zs=p->zs;
  double rs=p->rs;
  double g;
  double cosPhi, sinPhi;

  double xPlus, yPlus;
  double xMinus, yMinus;

  double dxPlus_dPhi, dxPlus_drs;
  double dyPlus_dPhi, dyPlus_drs;

  double dxMinus_drs, dxMinus_dPhi;
  double dyMinus_drs, dyMinus_dPhi;

  double xs, ys;		/* xs=rs*cosPhi, ys=rs*sinPhi */
  double f;			/* f=sqrt(1.0+4.0/g) */
  double h;			/* h=1/(g*g*sqrt(1+4/g)) */

  cosPhi = cos(phi); sinPhi = sin(phi);
  g = rs*rs+zs*zs+2.0*rs*zs*cosPhi;
  f = sqrt(1.0+4.0/g);
  h = 1.0/(g*g*f);

  xs = rs*cosPhi; ys=rs*sinPhi;

  xPlus = 0.5 * (zs+xs)* (1.0+f);
  yPlus = 0.5 *     ys * (1.0+f);

  xMinus = 0.5 *(zs+xs)* (1.0-f);
  yMinus = 0.5 *    ys * (1.0-f);

  dxPlus_drs  = cosPhi/2.0*(1.0+f) - 2.0*(zs+xs)*(rs+zs*cosPhi) * h;
  dxPlus_dPhi =-ys/2.0*    (1.0+f) + 2.0*(zs+xs)* zs*ys         * h;
  dyPlus_drs  = sinPhi/2.0*(1.0+f) - 2.0*    ys *(rs+zs*cosPhi) * h;
  dyPlus_dPhi = xs/2.0*    (1.0+f) + 2.0*    ys * zs*ys         * h;
  
  dxMinus_drs = cosPhi/2.0*(1.0-f) + 2.0*(zs+xs)*(rs+zs*cosPhi) * h;
  dxMinus_dPhi=-ys/2.0*    (1.0-f) - 2.0*(zs+xs)*zs*ys          * h;
  dyMinus_drs = sinPhi/2.0*(1.0-f) + 2.0*    ys *(rs+zs*cosPhi) * h;
  dyMinus_dPhi= xs/2.0*    (1.0-f) - 2.0*    ys *zs*ys          * h;

  return( - (-dxPlus_drs * dyPlus_dPhi + dyPlus_drs * dxPlus_dPhi)
	  + (-dxMinus_drs * dyMinus_dPhi + dyMinus_drs * dxMinus_dPhi) );

}

double gphmag(double zs, src_cld* src, double& Asecond)
{

  //calculate the gph magnification
  double A[13];

  double A00,A2rs2, A4rs4;
  double Ahplus, Afplus, Afcross;
  double Afirst,Agph;
  double zsi;
  Ahplus=Afplus=Afcross=0.0;

  for(int i=0;i<13;i++)
    {
      zsi = abs(zs + src->rs*cd(gphOffX[i],gphOffY[i]));
      A[i]=pacAmp(zsi);
    }

  A00=A[0];

  for(int i=1;i<5;i++) Ahplus+=A[i];
  Ahplus=0.25*Ahplus-A00;
  for(int i=5;i<9;i++) Afplus+=A[i];
  Afplus=0.25*Afplus-A00;
  for(int i=9;i<13;i++) Afcross+=A[i];
  Afcross=0.25*Afcross-A00;

  A2rs2=(16.0*Ahplus-Afplus)/3.0;
  A4rs4=0.5*(Afplus+Afcross)-A2rs2;

  Afirst=src->twoCoeff*A2rs2;
  Asecond=src->fourCoeff*A4rs4;
  Agph=A00+Afirst+Asecond;

  return Agph;
}


double leeInner(double x, void * params)
{
  struct innerParams * p = (struct innerParams *) params;
  src_cld* src = p->src;
  double zs = p->zs;
  double costheta = p->costheta;
  double x2=x*x;

  return (x2+2.0)/qAdd(x,2.0) * (1.0 - src->ldGamma * (1.0 - 1.5*sqrt(1.0 - (x2-2.0*x*zs*costheta+zs*zs)/sqr(src->rs))));
}

double leeOuter(double x, void * params)
{
  struct outerParams * p = (struct outerParams *) params;
  src_cld* src = p->src;
  double zs = p->zs;
  double u1,u2;
  double costheta=cos(x);
  double qud = sqrt(sqr(src->rs) - sqr(zs*sin(x)));

  struct innerParams ip={src,zs,costheta};

  if(zs<=src->rs)
    {
      u1=0.0;
      u2=zs*costheta + qud;
    }
  else if(x <= asin(src->rs/zs))
    {
      u1=zs*costheta - qud;
      u2=zs*costheta + qud;
    }
  else
    {
      u1=0.0;
      u2=0.0;
      return 0.0;
    }

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  
  double epsAbsolute=1.0e-5;
  double epsRelative=0.0;//1.0e-5;
  double result, error;

  gsl_function F;
  int status;

  gsl_set_error_handler_off();
  /* integrate to obtain magnification */
  F.function = &leeInner;
  F.params = &ip;

  status = gsl_integration_qags (&F, u1, u2, epsAbsolute, epsRelative, 10000, w, &result, &error);    
  if(status) {cerr << "wittFSPL.cpp: fsplMagnification: Error caught: " << status << endl; exit(1);}

  gsl_integration_workspace_free(w);

  return result;
}

double leeFSPLmagnification(double zs, src_cld* src)
{
  struct outerParams op={src,zs};

  double Agph, A2, mu;

  double epsAbsolute=1.0e-5;
  double epsRelative=0.0;//1.0e-5;
  double result, error;

  gsl_function F;
  int status;

  Agph = gphmag(zs,src,A2);

  if(abs(A2)<epsAbsolute)
    {
      mu=Agph;
    }
  else
    {

      gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  
      gsl_set_error_handler_off();
      /* integrate to obtain magnification */
      F.function = &leeOuter;
      F.params = &op;

      status = gsl_integration_qags (&F, 0.0, pi, epsAbsolute, epsRelative, 10000, w, &result, &error);    
      if(status) {cerr << "wittFSPL.cpp: fsplMagnification: Error caught: " << status << endl; exit(1);}

      gsl_integration_workspace_free(w);

      mu=2.0*result/(pi*sqr(src->rs));
    }

  return mu;
}
