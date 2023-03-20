#include <vector>
#include <string>
#include "structures.h"
#include <gsl/gsl_matrix.h>

using namespace std;

vector<double> mderrors(gsl_matrix* cov, int npar, int obs1, int obs2, double alpha, int ffp, struct event* Event, struct slcat *Sources, struct slcat *Lenses);
vector<double> piEdirerrors(gsl_matrix* cov, int npar, struct event* Event, struct slcat *Sources, struct slcat *Lenses);
