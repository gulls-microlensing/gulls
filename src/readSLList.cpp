#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<time.h>
#include<unordered_map>

#include "readSLList.h"
#include "split.h"
#include "columnCodes.h"

/*int data_index(string key)
{ // write a function to return the index of an element key in datakey
  
  for(int i=0;i<int(datakey.size());i++)
    {
      if (datakey[i] == key)
	{
	  return i;
	}
      else
	{
	  // Throw an error if the key is not found
	  throw invalid_argument("Key not found in datakey");
	}
    }
}
*/


int readSLList(int choosefield, bool src, struct filekeywords *Paramfile, struct slcat *sl)
{
  //int ncols = Paramfile->Nfilters + NDATAFIELDS;
  int nlist=0;
  int nlines=-1;
  int ncols;
  bool unlog_radius=false;
  bool radius_set=false;

  string line;
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");
  vector<string> sdata;
  vector<double> data;

  timespec timestart, timeend;
  clock_gettime(CLOCK_REALTIME,&timestart);

  ifstream sll;

  if(src) sll.open(Paramfile->sourcelist);
  else sll.open(Paramfile->lenslist);
  if(!sll)
    {
      cerr << __FUNCTION__ << "Error: Could not open " << (src?string("source"):string("lens")) << " list (" << (src?Paramfile->sourcelist:Paramfile->lenslist) << ")" << endl;
      return 0;
    }

  while(!sll.eof())
    {
      getline(sll,line);
      line = line.substr(0,line.find_first_of("#"));
      
      split(line,sdata);

      if(int(sdata.size())>=6)
	{
	  nlines++;
	  
	  int field;
	  double l,b;
	  double dl,db; //side lengths
	  string name;
	  
	  //read in the field information	  
	  field = stoi(sdata[0]);
	  
	  if(field != Paramfile->choosefield) continue;
	  
	  l = atof(sdata[1].c_str());
	  b = atof(sdata[2].c_str());
	  dl = atof(sdata[3].c_str());
	  db = atof(sdata[4].c_str());
	  name = sdata.back();

	  if(src)
	    {
	      Paramfile->avgl = l;
	      Paramfile->avgb = b;
	    }
	  
	  sl->field = field;
	  sl->start = nlist;
	  sl->l = l;
	  sl->b = b;
	  sl->dl = dl;
	  sl->db = db;
	  
	  //now read in the stars
	  string strfldname; 
	  if(src) strfldname = string(Paramfile->sourcedir) + name;
	  else strfldname = string(Paramfile->lensdir) + name;
	  
	  int nstarfieldlines=0;

	  if(choosefield<0 || choosefield==field)
	    {
	      cerr << "Loading field " << field << endl;
	      ifstream strfld(strfldname.c_str());
	      if(!strfld)
		{
		  cerr << __FUNCTION__ << ": Error: Could not open starfield " << field << " (" << strfldname << ")" << endl;
		  return 0;
		}
	      
	      while(!strfld.eof())
		{
		  getline(strfld,line);
		  
		  split(line,sdata);		  
		  
		  if(nstarfieldlines==0 && line.find_first_of(string("abcdghijklmnopqrstuvwxyzABCDGHIJKLMNOPQRSTUVWXYZ"))!=line.npos)
		    {
		      if(Paramfile->verbosity>=2)
			{
			  cout << "Reading SL header" << endl;
			}

		      ncols=sdata.size();
		      
		      //read in the magnitudes
		      if((src?Paramfile->sourcecolours:Paramfile->lenscolours))
			{
			  //magnitudes expressed as 1 mag + N colours
			  sl->magkey.push_back(sdata[0]);
			  sl->magdict[sdata[0]] = sl->magkey.size()-1;
			  for(int i=1;i<Paramfile->Nfilters;i++)
			    {
			      sl->magkey.push_back(sdata[i]);
			      sl->magdict[sdata[i]] = sl->magkey.size()-1;
			    }
			}
		      else
			{
			  //magnitudes all expressed as such
			  for(int i=0;i<Paramfile->Nfilters;i++)
			    {
			      sl->magkey.push_back(sdata[i]);
			      sl->magdict[sdata[i]] = sl->magkey.size()-1;
			    }
			}

		      //read in the other data
		      for(int i=Paramfile->Nfilters;i<ncols;i++)
			{
			  sl->datakey.push_back(sdata[i]);
			  sl->datadict[sdata[i]] = sl->datakey.size()-1;

			  //Adjust the default parameter numbers
			  if(sdata[i]==string("mul") || sdata[i]==string("mu_l")) sl->MUL = sl->datadict[sdata[i]];
			  if(sdata[i]==string("mub") || sdata[i]==string("mu_b")) sl->MUB = sl->datadict[sdata[i]];
			  if(sdata[i]==string("Vr") || sdata[i]==string("vr_bc")) sl->VR = sl->datadict[sdata[i]];
			  if(sdata[i]==string("UU") || sdata[i]==string("u") || sdata[i]==string("U")) sl->UU = sl->datadict[sdata[i]];
			  if(sdata[i]==string("VV") || sdata[i]==string("v") || sdata[i]==string("V")) sl->VV = sl->datadict[sdata[i]];
			  if(sdata[i]==string("WW") || sdata[i]==string("w") || sdata[i]==string("W")) sl->WW = sl->datadict[sdata[i]];
			  //if(sdata[i]==string("CL") || sdata[i]==string("phase")) CL = sl->datadict[sdata[i]];
			  //if(sdata[i]==string("TYP") || sdata[i]==string("phase")) CL = sl->datadict[sdata[i]];
			  //Currently user responsibility to check for logged values based on the key, may change later
			  if(sdata[i]==string("TEFF") || sdata[i]==string("Teff") || sdata[i]==string("log_Teff") || sdata[i]==string("logTeff")) sl->TEFF = sl->datadict[sdata[i]];
			  if(sdata[i]==string("LOGG") || sdata[i]==string("logg") || sdata[i]==string("log_g")) sl->LOGG = sl->datadict[sdata[i]];
			  if(sdata[i]==string("MASS") || sdata[i]==string("Mass") || sdata[i]==string("mass") || sdata[i]==string("final_mass")) sl->MASS = sl->datadict[sdata[i]];
			  if(sdata[i]==string("MBOL") || sdata[i]==string("logL") || sdata[i]==string("log_L")) sl->MBOL = sl->datadict[sdata[i]];

			  //This should mean that a radius column takes precedence over a logged column if both appear in the output
			  if(sdata[i]==string("RADIUS") || sdata[i]==string("Radius") || sdata[i]==string("radius"))
			    {
			      radius_set = true;
			      unlog_radius = false;
			      sl->RADIUS = sl->datadict[sdata[i]];
			    }
			  if(!radius_set && (sdata[i]==string("logRADIUS") || sdata[i]==string("log_R") || sdata[i]==string("log_radius")))
			    {
			      unlog_radius = true;
			      sl->RADIUS = sl->datadict[sdata[i]];
			    }

			  if(sdata[i]==string("LL") || sdata[i]==string("l") || sdata[i]==string("l(deg)")) sl->LL = sl->datadict[sdata[i]];
			  if(sdata[i]==string("BB") || sdata[i]==string("b") || sdata[i]==string("b(deg)")) sl->BB = sl->datadict[sdata[i]];
			  if(sdata[i]==string("RA") || sdata[i]==string("ra") || sdata[i]==string("RA2000.0")) sl->RA = sl->datadict[sdata[i]];
			  if(sdata[i]==string("DEC") || sdata[i]==string("dec") || sdata[i]==string("DEC2000.0")) sl->DEC = sl->datadict[sdata[i]];
			  if(sdata[i]==string("DIST") || sdata[i]==string("Dist") || sdata[i]==string("dist") || sdata[i]==string("distance")) sl->DIST = sl->datadict[sdata[i]];
			  if(sdata[i]==string("XX") || sdata[i]==string("x") || sdata[i]==string("X")) sl->XX = sl->datadict[sdata[i]];
			  if(sdata[i]==string("YY") || sdata[i]==string("y") || sdata[i]==string("Y")) sl->YY = sl->datadict[sdata[i]];
			  if(sdata[i]==string("ZZ") || sdata[i]==string("z") || sdata[i]==string("Z")) sl->ZZ = sl->datadict[sdata[i]];

			  
			}
		    }

		  if(Paramfile->verbosity>=2 && nstarfieldlines==0)
		    {
		      cout << "SL mag keys:" << endl;
		      for(int i=0;i<Paramfile->Nfilters;i++)
			{
			  cout << i << " " << sl->magkey[i] << endl;
			}
		      cout << "SL data keys:" << endl;
		      for(int i=0;i<sl->datakey.size();i++)
			{
			  cout << i << " " << sl->datakey[i] << endl;
			}		  
		    }
				  
		  //only read data lines
		  if(line.find_first_of(ignore)==line.npos)
		    {
		      split(line,data);
		      
		      if(int(data.size()) == ncols) 
			{
			  //read in the magnitudes
			  sl->mags.push_back(vector<double>(Paramfile->Nfilters));
			  if((src?Paramfile->sourcecolours:Paramfile->lenscolours))
			    {
			      //magnitudes expressed as 1 mag + N colours
			      sl->mags[nlist][0] = data[0];
			      for(int i=1;i<Paramfile->Nfilters;i++)
				{
				  sl->mags[nlist][i] = data[i] + data[0];
				}
			    }
			  else
			    {
			      //magnitudes all expressed as such
			      for(int i=0;i<Paramfile->Nfilters;i++)
				{
				  sl->mags[nlist][i] = data[i];
				}
			    }

			  //read in the other data
			  sl->data.push_back(vector<double>(sl->datakey.size()));
			  for(int i=Paramfile->Nfilters;i<ncols;i++)
			    {
			      sl->data[nlist][i-Paramfile->Nfilters] 
				= data[i];		      
			    }

			  nlist++;
			  nstarfieldlines++;
			}//end if(int(data.size()) == ncols) 
		    }//end if(line.find_first_of(ignore)==line.npos)
		}//end while(!strfld.eof())

	      sl->end = nlist-1;

	      strfld.close();

	    } //end choosefield

	}
      
    }

  sll.close();

  clock_gettime(CLOCK_REALTIME,&timeend);

  long nnanosec = timeend.tv_nsec - timestart.tv_nsec;
  double nsec = double((timeend.tv_sec - timestart.tv_sec) - (nnanosec<0?1:0))
    + double(nnanosec<0?nnanosec+1000000000:nnanosec)*1.0e-9;

  cout << "Reading " << (src?"sources":"lenses") << " took ";
  printf("%f seconds\n",nsec);
  
  return 1;

}

int getValidFields(struct filekeywords* Paramfile, struct slcat *Sources, struct slcat *Lenses, vector<double>* sfielddata)
{
  int sfok, sok, lok;
  sfok = sok = lok = 0;

  //check that there is a starfield
  for(int level=0;level<NUM_STARFIELD_LEVELS;level++)
    {
      if((*sfielddata)[level] > 0) sfok=1;
    }

  //check that there's a source list
  if(Sources->start >=0 && Sources->end>Sources->start)
    sok = 1;

  //check that there's a lens list
  if(Lenses->start >=0 && Lenses->end>Lenses->start)
    lok = 1;

  //if all elements needed are there - we have a valid field
  if(sfok && sok && lok)
    Paramfile->validFields.push_back(0); //field);

  cout << "sfok, sok, lok " << sfok << " " << sok << " " << lok << endl;
  cout << "Sources start, end " << Sources->start << " " << Sources->end << endl;
  cout << "Lenses start, end " << Lenses->start << " " << Lenses->end << endl;

  return Paramfile->validFields.size();
}
