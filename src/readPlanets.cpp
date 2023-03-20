#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<cstdlib>

#include "readPlanets.h"
#include "split.h"

int readPlanets(struct filekeywords *Paramfile, vector<struct pcat> *planets, char* instance, int choosefield)
{
  int ncols = NPLANETINPUT;
  int nlist=0;

  string line;
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");
  vector<double> data;

  ifstream pf;

  string fname = string(Paramfile->planetdir) + string(Paramfile->planetroot);

  if(choosefield>=0)
    {
      char field[20];
      sprintf(field,"%d",choosefield);
      fname += string(field) + string(".");
    }

  fname += string(instance);

  pf.open(fname.c_str());
  if(!pf)
    {
      cerr << "ERROR READING PLANETS FILE: " << fname << endl;
      return 0;
    }

  //read in the planets
  while(!pf.eof())
    {
      getline(pf,line);
      //remove any comments
      if(line.find_first_of(ignore)==line.npos)
	{
	  split(line,data);

	  if(int(data.size())==ncols)
	    {
	      planets->push_back(pcat(&data,NPLANETINPUT+NPLANETDERIV));
	      nlist++;
	    }
	}
    }

  return nlist;
}

