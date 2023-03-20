#include<cmath>

#include "cd.h"
#include "binaryMag.h"
#include "zroots2.h"
#include "integerPowers.h"
#include "pm.h"

/*double magnification(cd zs, int& nimgs, cd imgs[5], int& flag, lensParameters* l)
{
  const double eps=1.0e-7;
  cd z,df,dz;
  double J;
  double magi[5];
  cd solns[6];

  double magTotal=0.0;
  int nroots=0;

  cd dzsdzb;
  bool thesame;
  cd rem,swap;

  cd zsb = conj(zs);

  cd lc[6];

  flag=0;

  if(abs(zs)>4.0||nimgs<3||nimgs>5)
    {
      nimgs=3;
      imgs[0] = zs - zsb/((zsb-l->z1)*(zsb-l->z2));
      imgs[1] = l->z1 - l->m1/(zsb-l->z1-l->m2/(l->z2-l->z1));
      imgs[2] = l->z2 - l->m2/(zsb-l->z2-l->m1/(l->z1-l->z2));
    }

  //for each image
  for(int i=0;i<nimgs;i++)
    {
      z=imgs[i];

      //solve lens equation via newton raphson
      for(int j=0;j<100;j++)
	{
	  lensEquation(zs, z, df, J, dz, l);
	  z+=dz;

	  if(abs(dz)<eps)
	    {
	      solns[nroots]=z;
	      magi[nroots]=J;
	      nroots++;
	      break;
	    }
	}
    }

  //now check the solutions

  if(nroots>0)
    {
      //we have some solutions, check for duplicates

      imgs[0]=solns[0];
      magTotal=abs(magi[0]);
      nimgs=1;

      for(int i=1;i<nroots;i++)
	{
	  thesame=false;
	  for(int j=0;j<i;j++)
	    {
	      if(abs(solns[i]-solns[j])<2.0*eps) thesame=true;
	    }

	  if(!thesame)
	    {
	      imgs[nimgs]=solns[i];
	      magi[nimgs]=magi[i];
	      magTotal += abs(magi[i]);
	      nimgs++;
	    }
	}
    }

  //check if there are any other solutions
  if(nimgs<5)
    {
      coefficients(zs,l,lc); //get the coefficients of the 5th degree polynomial

      //deflate it
      for(int j=0;j<nimgs;j++)
	{
	  rem=lc[5-j];
	  for(int i=5-j-1;i>=0;i--)
	    {
	      swap = lc[i];
	      lc[i]=rem;
	      rem=swap+rem*imgs[j];
	    }
	}

      //now solve for the remaining images

      //if the remaining polynomial is quadratic, solve analytically
      if(nimgs==3)
	  qroot(lc[2], lc[1], lc[0], solns[1], solns[2]);
      else
	zroots(lc,5-nimgs,solns,1,"staticLens::nmagnification");

      //check if the new solutions are true solutions

      for(int i=1;i<=5-nimgs;i++)
	{
	  z=solns[i];
	  lensEquation(zs,z,df,J,dz,l);

	  //cout << abs(df) << " " << abs(J) << endl;

	  if(abs(df)<eps)
	    {
	      magi[nimgs] = J;
	      imgs[nimgs]=solns[i];
	      magTotal+=abs(J);
	      nimgs++;
	    }
	}
    }
  
  //check that we now have the right number of images
  if(!(nimgs==3 || nimgs==5))
    {
      flag=1;
    }

  //last resort if errors - solve 5th order with zroots
  if(flag)
    {
      flag = 0;
      nimgs=0;
      magTotal=0.0;
      coefficients(zs,l,lc);
      zroots(lc,5,solns,1,string("staticLens::nmagnification"));

      //check if the new solutions are true solutions

      for(int i=1;i<=5;i++)
	{
	  z=solns[i];
	  lensEquation(zs,z,df,J,dz,l);

	  if(abs(df)<eps)
	    {
	      magi[nimgs] = J;
	      imgs[nimgs]=solns[i];
	      magTotal+=abs(J);
	      nimgs++;
	    }
	}

      //check that we now have the right number of images
      if(!(nimgs==3 || nimgs==5))
	{
	  cerr << "Wrong number of images (" << nimgs << ") in nmagnification\n";
	  flag=1;
	}
    }

  //cout << nimgs << " " << magTotal << endl;
  return magTotal;
  }*/

double magnification(cd zs, int& nimgs, cd imgs[5], int& flag, lensParameters* l)
{
  const double eps=1.0e-7;
  cd z,df,dz;
  double J;
  double magi[5];
  cd solns[6];

  double magTotal=0.0;
  int nroots=0;

  cd dzsdzb;
  bool thesame;
  cd rem,swap;

  cd zsb = conj(zs);

  cd lc[6];

  flag = 0;
  nimgs=0;
  magTotal=0.0;
  coefficients(zs,l,lc);
  zroots(lc,5,solns,1,string("binaryMag::magnification"));

  //check if the new solutions are true solutions

  for(int i=1;i<=5;i++)
    {
      z=solns[i];
      lensEquation(zs,z,df,J,dz,l);

      if(abs(df)<eps)
	{
	  magi[nimgs] = J;
	  imgs[nimgs]=solns[i];
	  magTotal+=abs(J);
	  nimgs++;
	}
    }

  //check that we now have the right number of images
  if(!(nimgs==3 || nimgs==5))
    {
      //cerr << "Wrong number of images (" << nimgs << ") in nmagnification\n";
      flag=1;
    }

  return magTotal;

}

void coefficients(cd zs, lensParameters* l, cd coeffs[6])
{
  //lens equation coefficients in the frame where m1z2+m2z1=0
      
  cd g[6]={(l->z1*l->z2),(-l->z1-l->z2),1.0,0.0,0.0,0.0};
  cd zsz[6]={zs,(-1.0),0.0,0.0,0.0,0.0};
  cd zsb1=conj(zs)-l->z1;
  cd zsb2=conj(zs)-l->z2;
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
  poly5ConstMult(l->m1,part2,5);
      
  //part3
  poly5Mult(f1,g,part3,5);
  poly5ConstMult(l->m2,part3,5);
      
  //final coefficients
  poly5Add(coeffs,part2,5);
  poly5Add(coeffs,part3,5); 
      
}

void lensEquation(cd zs, cd& z, cd& df, double& jac, cd& dz,lensParameters* l)
{
  cd zb=conj(z);
  cd w, w1, w2;
  cd dfdzb;
  cd u;

  w1 = 1.0/(zb-l->z1);
  w2 = 1.0/(zb-l->z2);
  w = 1.0/zb;

  u = zb*w1*w2;

  df = z - zs - u;
  dfdzb = u*(w-w1-w2);

  jac = 1.0/(1.0 - sqr(real(dfdzb)) - sqr(imag(dfdzb)));

  dz = -(df + dfdzb * conj(df)) * jac;
}

bool analyticCaustic(double phi, lensParameters* l, cd cc[4], cd ca[4])
{

  //returns the position of the critical curve and caustic for a given phi
  //uses a coorditate system centred on the primary mass, but returns in a 
  //coordinate system at the centre of mass
  //cc and ca should be of size 4
  
  bool flag=false;

  double d=l->z2-l->z1;
  
  const double ot = 1.0/3.0;
  
  cd Wm, Wp, S, T, U;
  cd a_,b,c,c2,P,Q,R,V;
  double d_,d2;
  cd emiphi;
  cd Rm2,PpPoR,ToS;
  
  double m1=l->m1;
  double z1=l->z1;
  
  emiphi = cd(cos(-phi),sin(-phi));
  
  //coefficients
  
  //a_=-m1*x_*x_*emiphi;
  b=m1*d*emiphi;   
  a_=-d*b;
  //cout << "a_ = " << a_ << endl;
  b+=b; //times two - addition for efficincy b = 2*m1*x_*emiphi;
  //cout << "b = " << b << endl;
  c=-emiphi+d*d;
  //cout << "c = " << c << endl;
  c2 = c*c;
  d_=-d-d; //=-2*x_;
  //cout << "d = " << d << endl;
  d2=d_*d_;
  
  //first level auxilliary
  
  P = 12.0*a_+c2-3.0*d_*b;
  //cout << "P = " << P << endl;
  //cout << "P terms:  = 12.0*a_=" << 12.0*a_ << "+c2=" << c2 << "-3.0*d_*b" << -3.0*d_*b << endl;
  Q = 27.0*b*b + c*(-72.0*a_ + c2+c2) + 9.0*d_*(-b*c + 3.0*d_*a_);
  //cout << "Q = " << Q << endl;
  
  //if(!finite(real(P))||!finite(imag(P))) cout << "Bad P = " << P << ", ";
  //if(!finite(real(Q))||!finite(imag(Q))) cout << "Bad Q = " << Q << ", ";;

  //level two auxilliary
  
  if(abs(P)/abs(Q)<1.0e-3 && imag(P)==0 && imag(Q)==0)
    {
      V=abs(Q)*sqrt(1.0-4.0*P*sqr(P/Q));
      if(real(V*Q)<0.0) R=pow(P*sqr(P/Q),ot);
      else
	{
	  V = sqrt(-4.0*P*P*P + Q*Q);
	  R = pow(0.5*(Q+V),ot);
	}
    } 
  else
    {
      V = sqrt(-4.0*P*P*P + Q*Q);
      R = pow(0.5*(Q+V),ot);
    }
  
  //if(!finite(real(V))||!finite(imag(V))) cout << "Bad V = " << V << ", ";
  
  Rm2 = 1.0/(R*R);
  PpPoR = P/R + R;
  
  //if(!finite(real(R))||!finite(imag(R))) cout << "Bad R = " << R << ", ";
  //if(!finite(real(Rm2))||!finite(imag(Rm2))) cout << endl << "a_,b,c,d_,P,Q,V,R " << a_ << " " <<  b << " " <<  c << " " <<  d_ << " " <<  P << " " <<  Q << " " <<  V << " " <<  R << endl;
  //if(!finite(real(Rm2))||!finite(imag(Rm2))) cout << "Bad Rm2; R= " << R << " Rm2 = " << Rm2 << ", ";
  //if(!finite(real(PpPoR))||!finite(imag(PpPoR))) cout << "Bad PpPoR = " << PpPoR << ", ";
  //level three auxilliary
  
  S = sqrt(ot*(-c-c+0.75*d2+PpPoR));
  //if(!finite(real(S))||!finite(imag(S))) cout << "Bad S = " << S << ", ";

  T = -b-b + d_*(c - 0.25*d2);
  ToS = T/S;
  
  //if(!finite(real(T))||!finite(imag(T))) cout << "Bad T = " << T << ", ";
  //if(!finite(real(ToS))||!finite(imag(ToS))) cout << "Bad ToS = " << ToS << ", ";
  
  U = ot*(-4.0*c+1.5*d2-PpPoR);
  //if(!finite(real(U))||!finite(imag(U))) cout << "Bad U = " << U << ", ";
  
  //level four auxilliary
  
  Wm = sqrt(U-ToS);
  Wp = sqrt(U+ToS);
  
  //if(!finite(real(Wm))||!finite(imag(Wm))) cout << "Bad Wm = " << Wm << ", ";
  //if(!finite(real(Wp))||!finite(imag(Wp))) cout << "Bad Wp = " << Wp << ", ";
  
  //critical curves
  
  cc[0] = 0.5*(d-S-Wm) + z1; //d_=-2x_ => -0.5*d_=x_
  cc[1] = 0.5*(d-S+Wm) + z1; //last z1 term converts back to the standard 
  cc[2] = 0.5*(d+S-Wp) + z1; //coordinates of this class
  cc[3] = 0.5*(d+S+Wp) + z1;
  
  //consistency check

  for(int k=0;k<4;k++)
    {
      //if(!finite(real(cc[k]))||!finite(imag(cc[k]))) cout << "Bad cc = " << cc[k] << ", ";
      //cd tempcd, tempcd2;
      //double tempd;
      ca[k] = lensEquation(cc[k],l);
      //lensEquation(cc[k],ca[k],tempcd,tempd,tempcd2,l);
      //if(!finite(real(ca[k]))||!finite(imag(ca[k]))) cout << "Bad ca = " << ca[k] << ", ";
      //if(!finite(real(cc[k]))||!finite(imag(cc[k]))) cout << endl;
      flag=flag||checkCriticalCurve(cc[k],l);
    }
  
  //if there is a problem try with a different method
  if(flag)
    {
      flag=lageurCaustic(phi,l,cc,ca);
      //if(flag) cout << "flag" << endl;
    }

  //if there's still a problem, admit failure
  /*  if(flag) 
      {
      cout << "Error calculating caustic position " << check*check << endl;
      cout << "phi,t,d,q,cr: " << phi << " " << t << " " << d << " " << q << " " << causticRegime() << endl;
      }*/
  
  return flag;
  
}

bool newtonCaustic(double phi, lensParameters * l, cd cc[4], cd ca[4])
{
  for(int j=0;j<4;j++)
    {
      for(int k=0;k<100;k++)
	{
	  cd f = l->m1/sqr(cc[j]-l->z1) + l->m2/sqr(cc[j]-l->z2) 
	    - polar(1.0,phi);
	  cd df = -2.0 * ( l->m1/cube(cc[j]-l->z1) 
			     + l->m2/cube(cc[j]-l->z2) );
	  cc[j] -= f/df;
	  if(abs(f)<1.0e-10) break;
	}
      ca[j] = lensEquation(cc[j],l);
    }
  return true;
}

/*bool newtonCusp(double lambda, lensParameters* l, cd cusp[6])
{
  for(int j=0;j<6;j++)
    {
      for(int k=0; k<100; k++)
	{
	  cd doub = l->m1/sqr(cusp[j]-l->z1) + l->m2/sqr(cusp[j]-l->z2);
	  cd trip = l->m1/cube(cusp[j]-l->z1) + l->m2/cube(cusp[j]-l->z2);
	  cd qud = l->m1/hCube(cusp[j]-l->z1) + l->m2/hCube(cusp[j]-l->z2);
	  
	  cd f = sqr(trip)-lambda*cube(doub);
	  cd df = 6*trip*(sqr(doub)-qud);
	  cusp[j] -= f/df;
	  if(abs(f)<1.0e-10) break;
	}
    }
  return true;
  }*/

int lageurCusp(double lambda, lensParameters * l, cd cusp[7])
{
  cd coeffs[7],roots[7];
  
  cd a[4], b[3];
  int flag=0;
  double check;

  a[3] = cd(1,0);
  a[2] = cd(-3*(l->m1*l->z2 + l->m2*l->z1));
  a[1] = cd(3*(l->m1*sqr(l->z2) + l->m2*sqr(l->z1)));
  a[0] = cd(-(l->m1*cube(l->z2)+l->m2*cube(l->z1)));
  b[2] = cd(1,0);
  b[1] = a[2]*(2.0/3.0);
  b[0] = a[1]/3.0;

  cd a2[7], b2[5], b3[7];
  poly5Mult(a,3,a,3,a2);
  poly5Mult(b,2,b,2,b2);
  poly5Mult(b,2,b2,4,b3);
  poly5ConstMult(cd(-lambda,0),b3,6);
  poly5Add(a2,6,b3,6,coeffs);
  
  zroots(coeffs, 6, roots, 1, "binaryMag::lageurCusp");

  for(int i=1;i<=6;i++)
    {
      cd doub = l->m1/sqr(cusp[i]-l->z1) + l->m2/sqr(cusp[i]-l->z2);
      cd trip = l->m1/cube(cusp[i]-l->z1) + l->m2/cube(cusp[i]-l->z2);
      double check = abs(sqr(trip)-lambda*cube(doub));

      if(check>1e-10)
	{
	  flag |= (1<<i);
	}
      cusp[i]=roots[i];
    }

  return flag;
}

bool lageurCaustic(double phi, lensParameters * l, cd cc[4], cd ca[4])
{
  cd coeffs[5],roots[5];
  cd emiphi=cd(cos(-phi),sin(-phi));
  bool flag=false;
  double check;

  double d=abs(l->z2-l->z1);

  double d2=d*d;

  coeffs[0]=cd(l->m1*d2,0.0);
  coeffs[1]=cd(-2.0*l->m1*d,0.0);
  coeffs[2]=1.0-d2*emiphi;
  coeffs[3]=2.0*d*emiphi;
  coeffs[4]=-emiphi;

  zroots(coeffs, 4, roots, 1, "binaryMag::lageurCaustic");

  for(int i=0;i<4;i++) 
    {
      cc[i]=roots[i+1]+l->z1; //+z1 converts to correct coordinates
      check=abs(l->m1/sqr(cc[i]-l->z1) + l->m2/sqr(cc[i]-l->z2));
      if(abs(check-1.0)>1.0e-6)
	{
	  flag=true;
	}
      //cd tempcd, tempcd2;
      //double tempd;
      ca[i] = lensEquation(cc[i],l);
      //lensEquation(cc[i],ca[i],tempcd,tempd,tempcd2,l);
    }

  return flag;
}
