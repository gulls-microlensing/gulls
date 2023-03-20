#ifndef RADIAL_BASIS_FUNCTION_HEADER

#include<vector>
#include<complex>
#include<string>

#include<gsl/gsl_vector.h>

#define RBF_DEFAULT_C -1.0

using namespace std;

//class for computing radial basis function interpolations and their derivatives and integrals
class rbf{

 private:

  double c; //The smoothing scale
  double c2; //The square of the smoothing scale
  int ndim; //the number of dimensions
  int ndata; // the number of datapoints used for interpolation
  gsl_vector* lambda; //the coefficients of the interpolation
  vector<double> x;
  vector<double> scale; //the scaling for each axis

  int constructed;

 public:

  rbf()
    {
      constructed=0;
    }; //default constructor

  void clear()
  {
    x.clear();
    scale.clear();
    if(constructed) gsl_vector_free(lambda);
    constructed=0;
  };

  //Constructor - calculates the coefficients
  void construct_interpolation(vector<double>* _x, vector<double>* _y, double _c=RBF_DEFAULT_C);
  
  //Constructor - loads coefficients
  void construct_interpolation(vector<double> coeff, double _c);

  //recommended form as includes additional error checking
  double interpolate(vector<double>* _x, vector<double>* _y, double _c=RBF_DEFAULT_C)
    {
      if(constructed)
	{
	  cerr << __FILE__ << ": " << __FUNCTION__ << ": Error: Trying to construct an interpolation that is already constructed. Use clear before reconstructing." << endl;
	  exit(1);
	}
 
      construct_interpolation(_x, _y, _c);
      return c;
    };
  double interpolate(double _x[], double _y[], int _ndata, int _ndim=1, double _c=RBF_DEFAULT_C);

  ~rbf()
    {
      if(constructed) gsl_vector_free(lambda);
      constructed=0;
    };

  double interpolate(string input, double _c=RBF_DEFAULT_C);

  //Evaluation


  double evaluate(vector<double>* xpoint);
  double evaluate(double xpoint[], int ndim=1);
  double evaluate(double xpoint)
  {
    vector<double> xp;
    xp.push_back(xpoint);
    return evaluate(&xp);
  }

  double derivative(vector<double>* xpoint, int dim);
  double derivative(double xpoint[], int dim, int ndim=1);
  double derivative(double xpoint)
  {
    vector<double> xp;
    xp.push_back(xpoint);
    return derivative(&xp,0);
  }

  complex<double> derivative_fdf(vector<double>* xpoint, int dim);
  complex<double> derivative_fdf(double xpoint[], int dim, int ndim=1);
  complex<double> derivative_fdf(double xpoint)
  {
    vector<double> xp;
    xp.push_back(xpoint);
    return derivative(&xp,0);
  }

  //Integration

  //recommended form
  double integrate(vector<double>* xmin, vector<double>* xmax);
  double integrate(double xmin[], double xmax[], int _ndim=1);


};

#define RADIAL_BASIS_FUNCTION_HEADER
#endif
