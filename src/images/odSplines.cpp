#include<fstream>
#include<iostream>
#include<cmath>
#include<string>

#include<gsl/gsl_spline.h>

#include "split.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=3)
    {
      cerr << "Usage: ./odSplines <inputTable> <output>" << endl;
      exit(1);
    }

  int nl=81;
  int nb=81;
  int nd=70;

  string tablename = string(argv[1]);
  string outname = string(argv[2]);

  ifstream tab(tablename.c_str());
  if(!tab) 
    {
      cerr << "Could not open table file " << tablename << endl;
      exit(1);
    }
  
  ofstream out(outname.c_str());
  if(!out) 
    {
      cerr << "Could not open output file " << outname << endl;
      exit(1);
    }

  double* od = new double[nl*nb*nd];
  double* l = new double[nl];
  double* b = new double[nb];
  double* d = new double[nd];

  string line;
  vector<double> data;

  int ll=0, bb=0, dd=0; 

  //read in the tabulated optical depths

  while(!tab.eof())
    {
      getline(tab,line);
      split(line,data);

      if(line.compare(0,1,"#") && int(data.size())==6)
	{
	  od[nd*nb*ll + nd*bb + dd] = data[3];
	  l[ll] = data[0];
	  b[bb] = data[1];
	  d[dd] = data[2];
	  
	  dd++;
	  if(dd==nd)
	    {
	      dd=0; bb++;
	      if(bb==nb)
		{
		  bb=0; ll++;
		  if(ll==nl)
		    {
		      break;
		    }
		}
	    }
	} 
    }

  tab.close();

  cout << "table loaded" << endl;

  //calculate splines along the l direction

  vector<double> lout, bout;
  double* iod = new double[nd];
  double* d2y = new double[nd];
  gsl_spline* spline;
  
  for(int lo=0;lo<13;lo++)
    {
      lout.push_back(0.25*lo - 0.4);
      bout.push_back(0.25*lo - 3.2);
    }

  double ilstep = 1.0/(l[1]-l[0]);
  double ibstep = 1.0/(b[1]-b[0]);

  for(int lo=0;lo<int(lout.size());lo++)
    {
      ll = int((lout[lo] - l[0])*ilstep);
      double dl2 = l[ll+1] - lout[lo];
      double dl1 = lout[lo] - l[ll];

      cout << "bracket: " << l[ll] << " " << lout[lo] << " " << l[ll+1] << endl;

      for(int bo=0;bo<int(bout.size());bo++)
	{
	  //find the enclosing points in the tabulated function
	  //assumes regularly spaced tabulations
	  bb = int((bout[bo] - b[0])*ibstep);
	  double db2 = b[bb+1] - bout[bo];
	  double db1 = bout[bo] - b[bb];

	  double iden = ilstep*ibstep;

	  for(dd=0;dd<nd;dd++)
	    { 
	      iod[dd] = iden*( od[nd*nb*ll + nd*bb + dd] * dl2*db2
			       + od[nd*nb*(ll+1) + nd*bb + dd] * dl1*db2
			       + od[nd*nb*ll + nd*(bb+1) + dd] * dl2*db1
			       + od[nd*nb*(ll+1) + nd*(bb+1) + dd] * dl1*db1);
	    }

	  //compute the splines
	  int sno = int(bout.size())*lo + bo;
	  spline = gsl_spline_alloc(gsl_interp_cspline, nd);
	  gsl_interp_accel *acc = gsl_interp_accel_alloc();

	  gsl_spline_init(spline, d, iod, nd);
	  
	  for(int dd=0;dd<nd;dd++)
	    {
	      double d2y = gsl_spline_eval_deriv2(spline, d[dd], acc);
	      out << sno << " " << d[dd] << " " << iod[dd] << " " << d2y << "\n";
	    }
	  out << endl;

	  cout << "field " << sno  << " complete " << endl;
	  
	}
    }

  //clear up memory
  delete[] od;
  delete[] iod;
  delete[] d2y;
  delete[] l;
  delete[] b;
  delete[] d;
  
}
