#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include "split.h"
#include "constdefs.h"
//#include "astroFns.h"
#include <gsl/gsl_blas.h>
#include<cmath>

//#define PI 4.0*atan(1.0)
//#define TO_DEG  180.0 / (PI) 
//#define TO_RAD  PI / 180.0

//extern "C"
//{
  //void eq2gal(double in_coord1, double in_coord2, char dir, double *out_coord1, double *out_coord2);
//#include "astroFns.h"
//}

using namespace std;

void eq2gal(double in_coord1, double in_coord2, char dir, double *out_coord1, double *out_coord2){

/*  Input: (RA & DEC) || (GAL_LONG & GAL_LAT) in radians */
/*     Output: converted values in radians */
/*     dir = 'e' for Equatorial to Galactic; 'g' for Galactic to Equatorial */
/*  NOTE: eq2gal.m is WRONG - it has these the wrong way around */
 

 double a[] = {
-0.0548755604,+0.4941094279,-0.8676661490,
-0.8734370902,-0.4448296300,-0.1980763734,
-0.4838350155,+0.7469822445,+0.4559837762};

 double b[3];
 double ans[] = {0.0, 0.0, 0.0};

gsl_matrix_view Mat = gsl_matrix_view_array(a, 3, 3);

b[0] = cos(in_coord1)*cos(in_coord2);
b[1] = sin(in_coord1)*cos(in_coord2);
b[2] = sin(in_coord2);

gsl_vector_view Pos = gsl_vector_view_array(b, 3);

gsl_vector_view RPos = gsl_vector_view_array(ans, 3);


 /* printf("%f %f %f\n",b[0],b[1],b[2]); */

 switch(dir){
 case 'g':
   gsl_blas_dgemv(CblasNoTrans, 1.0, &Mat.matrix, &Pos.vector, 1.0, &RPos.vector);
   break;
 case 'e':
  gsl_blas_dgemv(CblasTrans, 1.0, &Mat.matrix, &Pos.vector, 1.0, &RPos.vector);
break;
default:
printf("unknown type of conversion\n");
break;
}

*out_coord1 = atan2(ans[1],ans[0]);
 *out_coord2 = atan(ans[2]/sqrt(ans[0]*ans[0] + ans[1]*ans[1]));
}

int main(int argc, char* argv[])
{

  if(argc!=3&&argc!=4) 
    {
      cerr << "Wrong number of arguments\n\nUsage:\n./sightLines <input> <output>\n" << endl;
      exit(1);
    }

  cout.precision(10);
  cout << PI << endl;

  string input, output;

  char dir='g';
  if(argc==4) dir='e';

  input = string(argv[1]);
  output = string(argv[2]);

  ifstream in(input.c_str());
  if(!in) 
    {
      cerr << "Could not open input file " << input << endl;
      exit(1);
    }

  ofstream out(output.c_str());
  if(!out) 
    {
      cerr << "Could not open output file " << output << endl;
      exit(1);
    }

  string line;
  vector<double> data;
  double ra,dec;

  //  cout << TO_DEG << endl;

  while(!in.eof())
    {
      getline(in,line);
      split(line,data);

      //      cout << data[0]*TO_RAD << " " << data[1]*TO_RAD << " " << ra << " " << dec << endl;
      ra=dec=0.0;

      eq2gal(data[0]*TO_RAD,data[1]*TO_RAD,dir,&ra,&dec);

      //      cout << ra << " " << dec << endl;

      ra = ra*TO_DEG;
      dec = dec*TO_DEG;

      out << ra << " " << dec << " " << data[0] << " " << data[1] << "\n";
    }

}
