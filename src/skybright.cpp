#include<iostream>

#include "astroFns.h"
#include "constants.h"

using namespace std;

int main()
{

  double jd = 2455378.726817;
  double ra = (17 * 51.0/60.0 + 53.1/3600.0)*15*pi/180.0;
  double dec = -(29 + 50.0/60.0 + 46.0/3600.0)*pi/180.0;
  double lat = -29.0146 *pi/180.0;
  double lon = -70.6926*pi/180.0;
  double height = 2380.0;
  double extcoeff = 0.0828;
  double zensky = 19.9;

  double deltaV, D, objmoondist, K;

  moon_sky_brightness(jd, ra , dec, lon, lat, height, extcoeff, zensky,  &deltaV, &D, &objmoondist, &K);

  cout << "Hello" << endl;
  cout << deltaV << " " << D << " " << objmoondist << " " << K << endl;
}
