#ifndef Z_R2
#define Z_R2

#include<cmath>
#include<cstdlib>
#include<complex>
#include<string>

using namespace std;
			
bool laguer(complex<double> a[], int m, complex<double> *x, int *its, string s);
double FMAX(double f1, double f2);

bool zroots(complex<double> a[], int m, complex<double> roots[], int polish, string s);
  

#endif /* z_R2 */
