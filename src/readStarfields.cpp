#include<vector>
#include<fstream>
#include<iostream>
#include<string>
#include<time.h>

#include "split.h"
#include "structures.h"
#include "definitions.h"

int readStarfields(int choosefield, struct filekeywords *Paramfile, vector<vector<vector<double> > >* sfield, vector<double>* sfd)
{

  //stored list of stellar magnitues
  //structure is sfield[besancon field number][level number][star number][band]

  static const int offset=0; //the first column (zero-indexed) with star mags

  timespec timestart, timeend;

  clock_gettime(CLOCK_REALTIME,&timestart);

  //int fdx;
  int nstars;
  string line;
  vector<string> data;
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");

  if(Paramfile->verbosity>1) cout << "Resize sfield" << endl;
  sfield->resize(NUM_STARFIELD_LEVELS);
  //sfield->push_back(vector<vector<double> >(NUM_STARFIELD_LEVELS));
  sfd->resize(NUM_STARFIELD_LEVELS, 0.0);

  //read in the list of starfields
  ifstream sfl(Paramfile->starfieldlist);
  if(!sfl)
    {
      cerr << __FUNCTION__ << "Error: Could not open starfield list (" << Paramfile->starfieldlist << ")" << endl;
      return 0;
    }

  if(Paramfile->verbosity>1) cout << "Start reading starfield list" << endl;
  while(!sfl.eof())
    {
      getline(sfl,line);
      split(line,data);

      if(int(data.size())>=4 && line.compare(0,1,"#"))
	{
	  int field, level;
	  double area;
	  string name;

	  //read in the field information	  
	  field = atoi(data[0].c_str());
	  level = atoi(data[1].c_str());
	  area = atof(data[2].c_str())*sqr(3600.0);
	  name = data.back();

	  //check for problems
	  if(level>=NUM_STARFIELD_LEVELS)
	    {
	      cerr << __FUNCTION__ << ": Error: level number (" << level << ") specified in fields list is larger than allowed (" << NUM_STARFIELD_LEVELS << ")" << endl;
	      return 0;
	    }
	  
	  if(Paramfile->verbosity>2) cout << "set sfd[level=" << level << "]" << " = " << area << endl; 
	  (*sfd)[level] = area;
	  nstars=0;

	  if(choosefield==field)
	    {
	      //now read in the stars
	      ifstream strfld((string(Paramfile->starfielddir) + name).c_str());
	      if(!strfld)
		{
		  cerr << __FUNCTION__ << ": Error: Could not open starfield (" << string(Paramfile->starfielddir) + name << ")" << endl;
		  return 0;
		}

	      if(Paramfile->verbosity>2) cout << "Read in stars from field " << field << endl;

	      while(!strfld.eof())
		{
		  getline(strfld,line);

		  if(line.find_first_of(ignore)==line.npos)
		    {
		      split(line,data);

		      if(int(data.size()) >= Paramfile->Nfilters+offset) 
			{
			  if(Paramfile->verbosity>3) cout << "Push back to (*sfield[level=" << level << "])" << endl;
			  if(Paramfile->verbosity>3) cout << "sfield->size()=" << sfield->size() << ", (*sfield)[level].size=" << (*sfield)[level].size() << endl;
			  (*sfield)[level].push_back(vector<double>(Paramfile->Nfilters));
			  for(int i=0;i<Paramfile->Nfilters;i++)
			    {
			      if(Paramfile->verbosity>3) cout << "Set (*sfield)[level][nstars=" << nstars << "][i=" << i << "], Nfilters=" << Paramfile->Nfilters << ", offset, data.size, i+offset=" << offset << " " << data.size() << " " << i+offset << endl;
			      (*sfield)[level][nstars][i] = 
				atof(data[i+offset].c_str());
			    } //end for i
			  nstars++;
			} //end if data size
		    } //end if ignore
		} //end while not eof

	      strfld.close();
	    } //end if choosefield
	}
      
    }

  sfl.close();

  clock_gettime(CLOCK_REALTIME,&timeend);

  long nnanosec = timeend.tv_nsec - timestart.tv_nsec;
  double nsec = double((timeend.tv_sec - timestart.tv_sec) - (nnanosec<0?1:0))
    + double(nnanosec<0?nnanosec+1000000000:nnanosec)*1.0e-9;

  cout << "Reading starfields took ";
  printf("%f seconds\n",nsec);

  return 1;
}
