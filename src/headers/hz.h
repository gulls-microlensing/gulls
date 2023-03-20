#ifndef KOPPARAPU_HZ

///////////////////////////////////////////////////////////////////////////////
//
//
//    C++ translation of Kopparapu et al's habitable zone estimates
//
//
///////////////////////////////////////////////////////////////////////////////

const int kopparapu_recent_venus = 0;
const int kopparapu_runaway_greenhouse = 1;
const int kopparapu_moist_greenhouse = 2;
const int kopparapu_maximum_greenhouse = 3;
const int kopparapu_early_mars = 4;
const int kopparapu_2AU_cloud_limit = 5;
const int kopparapu_1st_CO2_condensation_limit = 6;

//******************************************************************************
//This code calculates habitable zone 'fluxes' using the expression given in 
//the Kopparapu et al. paper. The corresponding output file is 'HZ_fluxes.dat'.
//It also generates a file 'HZ_coefficients.dat' that gives the coefficients for
//the analytical expression.
// Ravi kumar Kopparapu Feb 25 2012
//******************************************************************************

double hz_seff(double teff, int zone)
{

  //calculate the habitable zone by
  // a = ((L/Lsun)/Seff)^0.5 AU

  //Zones:
  //  0: Recent Venus
  //  1: Runaway greenhouse
  //  2: Moist greenhouse
  //  3: Maximum greenhouse
  //  4: Early Mars
  //  5: 2 AU cloud limit (Note: This limit is not given in the paper)
  //  6: 1st CO2 condensation limit (Note: This limit is not given in the paper)

  double seffsun[7]={1.7763, 1.0385, 1.0146, 0.3507, 0.3207, 0.2484, 0.5408};
  double a[7] = {1.4335e-4, 1.2456e-4, 8.1884e-5, 5.9578e-5, 5.4471e-5, 
		 4.2588e-5, 4.4499e-5};
  double b[7] = {3.3954e-9, 1.4612e-8, 1.9394e-9, 1.6707e-9, 1.5275e-9, 
		 1.1963e-9, 1.4065e-10};
  double c[7] = {-7.6364e-12, -7.6345e-12, -4.3618e-12, -3.0058e-12, 
		 -2.7481e-12, -2.1709e-12, -2.2750e-12};
  double d[7] = {-1.1950e-15, -1.7511E-15, -6.8260e-16, -5.1925e-16, 
		 -4.7474e-16, -3.8282e-16, -3.3509e-16};
  int i=zone;

  double tstar = teff - 5780.0; //delta t relative to the sun

  //Calculating HZ fluxes for stars with 2600 K < T_eff < 7200 K
  return seffsun[i] + tstar*(a[i] + tstar*(b[i] + tstar*(c[i] + tstar*d[i])));

}

#define KOPPARAPU_HZ
#endif
