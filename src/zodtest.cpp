#include<iostream>
#include<cstdlib>
#include<cmath>

#include "zodiacalLight.h"

using namespace std;

int main(int argc, char* argv[])
{

  double min=0.97;
  double max=2.0;

  if(argc>2)
    {
      min = atof(argv[1]);
      max = atof(argv[2]);
    }

  zodiacalLight z;

  //z.outputSpectrum();

  vector<double> lambda;
  vector<double> eff;
  
  for(double l=0.6;l<2.7;l+=0.01)
    {
      lambda.push_back(l);
      if(l>=min&&l<=max) eff.push_back(1.0);
      else eff.push_back(0.0);
    }

  z.set_bandpass(lambda,eff);

  for(double x=-20;x<=380;x+=1)
    {
      for(double y=-110;y<110;y+=1)
	{
	  cout << x << " " << y << " " << z.get_mag(x,y) << " " << z.get_mag(x,5) << endl;
	}
      cout << endl;
    }

  return 0;
}

