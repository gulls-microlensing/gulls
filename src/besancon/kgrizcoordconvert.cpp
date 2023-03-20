#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cstdlib>
#include<cmath>

#include "coords.h"
#include "split.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=2)
    {
      cerr << "Usage: ./coordtest <input>" << endl;
      exit(1);
    }

  coords c;
  int nfilters=15;
  int nfiltold=10;
  int dfilt=nfilters-nfiltold;

  ifstream in(argv[1]);
  if(!in)
    {
      cerr << "Could not open input file " << argv[1] << endl;
      exit(1);
    }

  while(!in.eof())
    {
      string line;
      getline(in,line);

      if(in.eof() && line.size()==0) break;

      if(line.find_first_of("#")!=string::npos)
	{
	  cout << line << endl;
	  continue;
	}

      vector<string> data;
      split(line,data);

      
      double l = atof(data[26+dfilt].c_str())*pi/180.0;
      double b = atof(data[27+dfilt].c_str())*pi/180.0;
      double D = atof(data[30+dfilt].c_str());
      double U = atof(data[13+dfilt].c_str());
      double V = atof(data[14+dfilt].c_str());
      double W = atof(data[15+dfilt].c_str());
      double X = atof(data[31+dfilt].c_str());

      //Correct for a bug in the Besancon model
      if(X>8.5) V = -(2*226.4)-V; 

      //printf("l=%g b=%g D=%g U=%g V=%g W=%g\n",l,b,D,U,V,W);

      vector<double> rvmuadlb;
      c.uvw2rvmuadlb(l,b,D,U,V,W,&rvmuadlb);
      //printf("%12.6f %12.6f %12.6f %12.6f %12.6f ",rvmuadlb[0],rvmuadlb[1],rvmuadlb[2],rvmuadlb[3],rvmuadlb[4]);

      char tmp[30];
      sprintf(tmp,"%7.2f",V); data[14+dfilt] = string(tmp);
      sprintf(tmp,"%7.3f",rvmuadlb[3]); data[10+dfilt] = string(tmp);
      sprintf(tmp,"%7.3f",rvmuadlb[4]); data[11+dfilt] = string(tmp);
      sprintf(tmp,"%7.2f",rvmuadlb[0]); data[12+dfilt] = string(tmp);

      for(int i=0;i<int(data.size());i++)
	{
	  if(i>0) cout << " ";
	  cout << data[i];
	}
      cout << "\n";
    }
      
}
