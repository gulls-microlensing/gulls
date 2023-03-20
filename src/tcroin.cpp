#include<cmath>
#include<iostream>
#include<cstdlib>

#include "croin.h"
#include "constants.h"

using namespace std;

int main()
{
  double ucmax, t0, u0, xc, yc;

  double s=1.7;
  double q=1e-5;
  double tc=4.2;
  double uc=0.3;
  double tE=0.5;
  double alpha = 74;
  double ca=cos(alpha*pi/180);
  double sa=sin(alpha*pi/180);

  usecroin(s, q, tc, tE, uc, alpha, &t0, &u0, &ucmax);
  cout << "#" << t0 << " " << u0 << " " << xc << " " << yc << endl;
  //exit(1);

  double x1,y1,x2,y2;

  for(double t=-1;t<1;t+=0.01)
    {
      double tau = (t-t0)/tE;
      double tauc = (t-tc)/tE;
      
      x1 = tau*ca - u0*sa;
      y1 = tau*sa + u0*ca;

      x2 = tauc*ca - uc*ucmax*sa + xc;
      y2 = tauc*sa + uc*ucmax*ca + yc;

      cout << t << " " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
    }
}
