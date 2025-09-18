#include "readParamfile.h"
#include<unordered_map>
#include<fstream>
#include<iostream>

#include "split.h"

using namespace std;

void readParamfile(string v_file, struct filekeywords *Paramfile){
  cout << "Parameter file: " << v_file << endl;
  //for reading and later combing keywords
  //S List of the keywords in the parameter files. Now requiring a path to the 'paths file' and the file name

  unordered_map<string,int> pfdefault;
  
  //The expected parameters and their defaults
  unordered_map<string,string> pfile = {
    {"OBSERVATORY_DIR",""},
    {"OBSERVATORY_LIST",""},
    {"SET_RANDOM_SEED_TO_CLOCK","1"},
    {"RANDOM_SEED","1"},
    {"SIMULATION_ZERO_TIME",""},
    {"WEATHER_PROFILE_DIR",""},     
    {"RUN_NAME",""},
    {"OUTPUT_DIR",""},
    {"FINAL_DIR",""},          //not used in cpp, but in postprocessing
    {"EXECUTABLE",""},          //not used in cpp, but used in launch scripts
    {"RATES_FILE",""},
    {"PRINCIPLE_OBSERVATORY","0"}, 
    {"OUTPUT_LC","0"},
    {"STARFIELD_DIR",""},
    {"STARFIELD_LIST",""},    
    {"SOURCE_DIR",""},
    {"SOURCE_LIST",""},
    {"SOURCE_COLOURS","0"},     
    {"LENS_DIR",""},
    {"LENS_LIST",""},
    {"LENS_COLOURS","0"},           
    {"PLANET_DIR",""},
    {"PLANET_ROOT",""},
    {"NFILTERS",""},
    {"AMIN",""},    
    {"LARGEPSFMAG",""},
    {"OUTPUT_IMAGES","0"},
    {"PRETTY_PICS","0"},      
    {"PRETTY_PICS_DIMENSIONS","256"},
    {"MIN_CHISQUARED",""},     
    {"OUTPUT_ONERR","0"},
    {"OUTPUT_ONDET","0"},
    {"OUTPUT_ONALL","0"},    
    {"PARALLAX","1"},
    {"LENS_LIGHT","1"},
    {"REPEAT_SEQUENCE","0"}, 
    {"OBS_GROUPS","(ALL)"},
    {"OBS_GROUP_NAMES",""},   //not used in cpp at the moment, but is in postprocessing
    {"NUM_SIM_DAYS","2010"},
    {"U0MAX","3"},           
    {"ERROR_SCALING","0"},   
    {"SUBRUNSIZE",""},
    {"LC_GEN","1"},
    {"LD_GAMMA","0.00"},     
    {"VBM_RELTOL","1.0e-6"},
    {"VBM_ABSTOL","1.0e-4"},
    {"LC_TIMEOUT","60.0"},
    {"MULTIPLE_SOURCES","0"},
    {"MULTIPLE_LENSES","0"}
  };

  //For testing which parameters are at their default values
  for(auto it=pfile.begin(); it!=pfile.end(); it++)
    {
      pfdefault[it->first] = 1;
    }
  
  int errsum=0; //use this to check that critical keywords are included in the parameter file
                //test read_config_var return to set defaults if the keyword is not

  //getenv here for the path to get base_path, tack on to the beginning of everything, should not be bad this way, we can ditch paths.txt or whatever, will make everything simpler.

  char const* tmp;

  Paramfile->basedir=string("");
  
  if(!(tmp = getenv("GULLS_BASE_DIR")))
    {
      cout << "GULLS_BASE_DIR environment variable not set" << endl;
      exit(1);
    }
  else
    {
      Paramfile->basedir = string(tmp);
    }
  cout << "GULLS_BASE_DIR:" << Paramfile->basedir << endl;
  
  if(tmp = getenv("GULLS_STARS_DIR"))
    {
      Paramfile->starsdir = string(tmp); 
    }
  else
    {
      cout << "GULLS_STARS_DIR environment variable not set, assuming it is the same as GULLS_BASE_DIR" << endl;
      Paramfile->starsdir = Paramfile->basedir;
    }

  //Read in all the parameters
  ifstream f;
  
  f.open(v_file.c_str());
  if(!f)
    {
      cout << __FUNCTION__ << "Error: Could not open parameter file (" << v_file << ")" << endl;
      exit(1);
    }

  string key, value;
  string line, lineorig;
  vector<string> data;
  
  while(!f.eof())
    {
      getline(f,lineorig);
      line = lineorig.substr(0,lineorig.find_first_of("#"));
      
      split(line,data,"=");

      if(data.size()<2)
	{
	  if(Paramfile->verbosity>2) cout << "Skipping line in parameter file: " << lineorig << endl;
	  continue;
	}

      if(data.size()>2)
	{
	  if(Paramfile->verbosity>2) cout << "Comment in parameter file: " << lineorig << endl;
	}

      
      key = data[0]; trim(key);
      value = data[1]; trim(value);

      if(Paramfile->verbosity>2) cout << "Read from parameter file: " << key << "=" << value << endl;
      
      if(pfile.find(key)==pfile.end())
	{
	  cout << "WARNING: Item " << key << "=" << value << " in the parameter file is not recognized, it will have no effect on the run" << endl;
	}
      pfile[key] = value;
      pfdefault[key] = 0;
    }

  //Check for missing parameters without defaults
  for(auto it = pfile.begin(); it!=pfile.end(); it++)
    {
      if(it->second.length()==0)
	{
	  cout << __FUNCTION__ << ": ERROR: Required parameter " << it->first << " is not set in the parameter file." << endl;
	  errsum++;
	}
    }

  if(errsum>0)
    {
      cout << __FUNCTION__ << ": There were " << errsum << " missing values in the parameter file (" << v_file << "). Exiting" << endl;
      exit(1);
    }
	  

  //All the data that we need should be in pfile, start assigning it to variables (and parse it further where needed)

  //Directories and other input files
  Paramfile->run_name = pfile["RUN_NAME"];
  Paramfile->outputdir = pfile["OUTPUT_DIR"] + Paramfile->run_name + string("/");

  Paramfile->obsdir = Paramfile->basedir + pfile["OBSERVATORY_DIR"];
  Paramfile->obslist = Paramfile->obsdir + pfile["OBSERVATORY_LIST"];
  Paramfile->weatherprofiledir = Paramfile->basedir + pfile["WEATHER_PROFILE_DIR"];

  Paramfile->starfielddir = Paramfile->starsdir + pfile["STARFIELD_DIR"];
  Paramfile->starfieldlist = Paramfile->starfielddir + pfile["STARFIELD_LIST"];

  Paramfile->sourcedir = Paramfile->starsdir + pfile["SOURCE_DIR"];
  Paramfile->sourcelist = Paramfile->sourcedir + pfile["SOURCE_LIST"];

  Paramfile->lensdir = Paramfile->starsdir + pfile["LENS_DIR"];
  Paramfile->lenslist = Paramfile->lensdir + pfile["LENS_LIST"];

  Paramfile->planetdir = Paramfile->basedir + pfile["PLANET_DIR"];
  Paramfile->planetroot = pfile["PLANET_ROOT"];

  
  //Basic setup
  Paramfile->setseedtoclock = stoi(pfile["SET_RANDOM_SEED_TO_CLOCK"]);
  if(Paramfile->setseedtoclock==0)
    Paramfile->Seed = stoi(pfile["RANDOM_SEED"]);
  Paramfile->simulation_zerotime = stod(pfile["SIMULATION_ZERO_TIME"]);
  

  //Additional
  if(pfile.find("PRINCIPAL_OBSERVATORY")!=pfile.end())
    Paramfile->principle_observatory = stoi(pfile["PRINCIPAL_OBSERVATORY"]);
  else
    Paramfile->principle_observatory = stoi(pfile["PRINCIPLE_OBSERVATORY"]);

  Paramfile->outputLightcurve=stod(pfile["OUTPUT_LC"]);
  Paramfile->sourcecolours=stoi(pfile["SOURCE_COLOURS"]);
  Paramfile->lenscolours=stoi(pfile["LENS_COLOURS"]);
  Paramfile->Nfilters = stoi(pfile["NFILTERS"]);
  Paramfile->Amin = stod(pfile["AMIN"]);
  Paramfile->large_psf_mag = stod(pfile["LARGEPSFMAG"]);
  Paramfile->outputImages = stoi(pfile["OUTPUT_IMAGES"]);
  Paramfile->prettypic = stoi(pfile["PRETTY_PICS"]);
  Paramfile->minChiSquared = stod(pfile["MIN_CHISQUARED"]);
  Paramfile->outputOnErr = stoi(pfile["OUTPUT_ONERR"]);
  Paramfile->outputOnDet = stoi(pfile["OUTPUT_ONDET"]);
  Paramfile->outputOnAll = stoi(pfile["OUTPUT_ONALL"]);
  Paramfile->pllxMultiplyer = stod(pfile["PARALLAX"]);
  Paramfile->lenslight = stoi(pfile["LENS_LIGHT"]);
  Paramfile->identicalSequence = stoi(pfile["REPEAT_SEQUENCE"]);
  Paramfile->NUM_SIM_DAYS = stoi(pfile["NUM_SIM_DAYS"]);
  Paramfile->u0max = stod(pfile["U0MAX"]);
  Paramfile->error_scaling = stoi(pfile["ERROR_SCALING"]);
  Paramfile->SUBRUNSIZE = stoi(pfile["SUBRUNSIZE"]);
  Paramfile->LC_GEN = stoi(pfile["LC_GEN"]);
  Paramfile->LD_GAMMA = stod(pfile["LD_GAMMA"]);
  Paramfile->vbm_reltol = stod(pfile["VBM_RELTOL"]);
  Paramfile->vbm_tol = stod(pfile["VBM_ABSTOL"]);  
  Paramfile->lc_timeout = stod(pfile["LC_TIMEOUT"]);
  Paramfile->multiple_sources = stoi(pfile["MULTIPLE_SOURCES"]);
  Paramfile->multiple_lenses = stoi(pfile["MULTIPLE_LENSES"]);
  
  //Obsgroups
  Paramfile->obsgroupstr = pfile["OBS_GROUPS"];

  //Pretty pic dimensions
  size_t pos;
  pos = pfile["PRETTY_PICS_DIMENSIONS"].find_first_of(",xX:");
  if(pos==string::npos)
    {
      Paramfile->prettypicDimX = stoi(pfile["PRETTY_PICS_DIMENSIONS"]);
      Paramfile->prettypicDimY = Paramfile->prettypicDimX;
    }
  else
    {
      Paramfile->prettypicDimX = stoi(pfile["PRETTY_PICS_DIMENSIONS"].substr(0,pos));
      Paramfile->prettypicDimY = stoi(pfile["PRETTY_PICS_DIMENSIONS"].substr(pos+1));
    }

  cout << "Input file " << v_file << " PARSED" << endl;
  cout << "----------------------------------" << endl;
  for(auto it = pfile.begin(); it!=pfile.end(); it++)
    {
      cout << it->first << "=" << it->second;
      if(pfdefault[it->first]==1) cout << " [DEFAULT]";
      cout << endl;
    }
  cout << "----------------------------------" << endl;
  
}
