#ifndef FISHERINVERSIONHEADER

#include<string>
#include<vector>
#include <gsl/gsl_matrix.h>

using namespace std;

void fisherInversion(vector<double> dF, vector<double> err, int npar, int ndata, string filename, string header, string values, gsl_matrix** covmatrix=NULL, gsl_matrix** covinverse=NULL);

#define FISHERINVERSIONHEADER
#endif
