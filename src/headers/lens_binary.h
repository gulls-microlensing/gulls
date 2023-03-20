#ifndef LENS_BINARY
#define LENS_BINARY

#include "integerPowers.h"
#include "zroots2.h"
#include "pm.h"
#include "lens_base.h"
#include "binaryMag.h"

using namespace std;

//lens specific constants
//static const int nImagesMax=5;

class lens_binary : public lens_base
{

 private:

  lensParameters ll;

 public: 

  double z1, z2;     //lens positions for the binary lens
  double m1, m2;     //relative lens masses

  lens_binary()
    {
      double d=0.5;
      double q=0.5;

      ll = lensParameters(d,q);

      m1=1.0/(1.0+q);
      m2=q*m1;
      z1=-d/(1.0+q);
      z2=-q*z1;

      setBinaryLensConfig();
      checkForChanges();

    };

  lens_binary(double d, double q)
    {

      m1=1.0/(1.0+q);
      m2=q*m1;
      z1=-d/(1.0+q);
      z2=-q*z1;

      ll = lensParameters(d,q);

      if(d<1.0e-8 || q<1.0e-10)
	{
	  cerr << "Warning, using single lens approximation" << endl;
	  setSingleLensConfig();
	}
      else
	{
	  setBinaryLensConfig();
	}

      checkForChanges();
    };

  lens_binary(double m1_, double z1_, double m2_, double z2_)
    {
      m1=m1_; m2=m2_; z1=z1_; z2=z2_;

      ll = lensParameters(m1,m2,z1,z2);

      if(m2_<1.0e-10 || abs(z2_-z1_)<1.0e-10)
	{
	  cerr << "Warning, using single lens approximation" << endl;
	  setSingleLensConfig();
	}
      else
	{
	  setBinaryLensConfig();
	}

      checkForChanges();
    };

  ~lens_binary(){};

  inline void lensEquation(cd image, cd& source)
  {
    source=image - m1/conj(image-z1) - m2/conj(image-z2);
  }

  inline cd lensEquation(cd image)
  {
    return image - m1/conj(image-z1) - m2/conj(image-z2);
  }

  inline void verboseLensEquation(cd& z, cd& zs, cd& df, double& jac)
  {
    cd zb=conj(z);
    cd w, w1, w2;
    cd dfdzb;
    cd u;
    double dr2, di2;
    
    w1 = 1.0/(zb-z1);
    w2 = 1.0/(zb-z2);
    w = 1.0/zb;
    
    u = zb*w1*w2;
    
    df = z - zs - u;
    dfdzb = u*(w-w1-w2);
    dr2 = sqr(real(dfdzb));
    di2 = sqr(imag(dfdzb));
    
    jac = 1.0/(1.0 - dr2 - di2);
    
    //dz = -(df + dfdzb * conj(df)) * jac;
  }

  inline void getleCoeffs(cd zs, cd* coeffs)
  {
    //lens equation coefficients in the frame where m1z2+m2z1=0
    
    cd zsb=conj(zs);
      
    cd g[6]={(z1*z2),(-z1-z2),1.0,0.0,0.0,0.0};
    cd zsz[6]={zs,(-1.0),0.0,0.0,0.0,0.0};
    //cd z_[6]={0.0,1.0,0.0,0.0,0.0,0.0};
    cd zsb1=zsb-z1;
    cd zsb2=zsb-z2;
    cd f1[6], f2[6], part1[6], part2[6], part3[6];
  
    //f1
    poly5ConstMult(zsb1,g,f1,5);
    f1[1]+=1.0;
  
    //f2
    poly5ConstMult(zsb2,g,f2,5);
    f2[1]+=1.0;
  
    //part1
    poly5Mult(f1,f2,part1,5);
    poly5Mult(part1,zsz,coeffs,5);
  
    //part2
    poly5Mult(f2,g,part2,5);
    poly5ConstMult(m1,part2,5);
  
    //part3
    poly5Mult(f1,g,part3,5);
    poly5ConstMult(m2,part3,5);
  
    //final coefficients
    poly5Add(coeffs,part2,5);
    poly5Add(coeffs,part3,5);
  }

  inline int solveLensEquation(cd zs, int& nimgs, cd* images)
  {
    int flag;

    if(nImagesMax==5)
      {
	magnification(zs,nimgs,images,flag,&ll);
      }
    else
      {
	//assume we are in the single lens approximation
	nimgs=2;
	zs = zs - cd(-m1*abs(z2-z1),0);
	images[0]=0.5*sqrt(zs*zs+cd(4,0));
	images[1]=0.5*zs-images[0];
	images[0]+=0.5*zs;
	flag=0;
      }

    return flag;
  }

  inline double ssolveLensEquation(cd zs, int& nimgs, cd* images, int* flag)
  {
    double mag;
    if(nImagesMax==5)
      {
	int tmpflag;
	mag = magnification(zs,nimgs,images,tmpflag,&ll);
	*flag = tmpflag;
      }
    else
      {
	nimgs=2;
	images[0]=0.5*sqrt(zs*zs+cd(4,0));
	images[1]=0.5*zs-images[0];
	images[0]+=0.5*zs;
	double u2=sqr(abs(zs));
	mag = (u2+2)/sqrt(u2*(u2+4));
	*flag=0;
      }
    return mag;
  }

  void setBinaryLensConfig()
  {
    nImagesMax=5;
    nleCoeffs=6;
    nPossibleImageConfigs=2;
    imageConfig.resize(nPossibleImageConfigs);
    imageConfig[0]=3; 
    imageConfig[1]=nImagesMax;

    shift=(z2-z1)*(m1-m2)/(m1+m2);
  }

  void setSingleLensConfig()
  {
    //set to the parameters of a single lens for cases where certain 
    //masses or separations are zero
    nImagesMax=2;
    nleCoeffs=6;  //the lens equation should still be calculated correctly
    nPossibleImageConfigs=1;
    imageConfig.resize(nPossibleImageConfigs);
    imageConfig[0]=nImagesMax; 

    shift=0.0;
  }

  void checkForChanges()
  {
    static double oldm1;
    static double oldm2;
    static double oldz1;
    static double oldz2;

    if(!(m1==oldm1 && m2==oldm2 && z1==oldz1 && z2==oldz2))
      {
	//if there have been changes make it known
	//cout << "changed " << m1-oldm1 << " " << m2-oldm2 << " " << z1-oldz1 << " " << z2-oldz2 << endl;
	changeFlag=true;
	oldm1=m1; oldm2=m2;
	oldz1=z1; oldz2=z2;
      }
    //Never change the change flag to false - some dependent processes
    //may not have learned of the change yet!
  }

  bool changed() 
  {
    //default
    return changeFlag;
  }


};


#endif /*LENS*/
