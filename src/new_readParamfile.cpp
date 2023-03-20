#include "readParamfile.h"
#include <string>
#include <stdio.h>
using namespace std;

void readParamfile(char *v_file, struct filekeywords *Paramfile){
  //p_file is path and location to a text file that contains all the paths to different files, e.g. where the paramter file is, where the weather files are, etc.
  char *p_file;
  char str[1000];
  char str0[5];
  char str1[5];
  char str2[5];
  char str3[15];
  char str4[5];
  char str5[15];
  char str6[15];
  char str7[30];
  char str8[5];
  char str9[5];
  char str10[5];
  char str11[30];
  char str12[20];
  char str13[5];
  char str14[5];
  char str15[5];
  char str16[5];
  char str17[20];
  char str18[5];
  char str19[20];
  char str20[20];
  //S List of the keywords in the parameter files. Now requiring a path to the 'paths file' and the file name
  const char *keywords[39] = {"OBSERVATORY_DIR", "OBSERVATORY_LIST",             //0,1
			      "SET_RANDOM_SEED_TO_CLOCK", "RANDOM_SEED",         //2,3
			      "SIMULATION_ZERO_TIME", "WEATHER_PROFILE_DIR",     //4,5
			      "RUN_NAME", "OUTPUT_DIR", "PRINCIPLE_OBSERVATORY", //6,7,8
			      "OUTPUT_LC", "STARFIELD_DIR", "STARFIELD_LIST",    //9,10,11
			      "SOURCE_DIR", "SOURCE_LIST", "SOURCE_COLOURS",     //12,13,14
			      "LENS_DIR", "LENS_LIST", "LENS_COLOURS",           //15,16,17 
			      "PLANET_DIR","PLANET_ROOT", "NFILTERS", "AMIN",    //18,19,20,21
			      "LARGEPSFMAG","OUTPUT_IMAGES", "PRETTY_PICS",      //22,23,24
			      "PRETTY_PICS_DIMENSIONS", "MIN_CHISQUARED",        //25,26
			      "OUTPUT_ONERR", "OUTPUT_ONDET", "OUTPUT_ONALL",    //27,28,29
			      "PARALLAX","LENS_LIGHT","REPEAT_SEQUENCE",         //30,31,32
			      "OBS_GROUPS","NUM_SIM_DAYS","U0MAX",               //33,34,35
			      "PATHS_DIR","PATHS_FILE","BASE_DIR"};              //36,37,38

  //keywrods in the path file, I believe all are duplicate except for base_dir, so let's only read in that one
  const char *path_keywords[7] = {"BASE_DIR","SCRIPT_DIR",        //0,1 
				  "PARAM_DIR","SRC_DIR",          //2,3
				  "PARAM_FILE_DIR","OBSERV_DIR",  //4,5
				  "WEATHER_DIR"};                 //6
  int errsum=0; //use this to check that critical keywords are included in the parameter file
                //test read_config_var return to set defaults if the keyword is not
 
  //read the parameters

  //starting with the paths file first, as it will be needed for later stuff
  //NOTE: this still requires the entire path for the paramfile as an arguement, for v_file I mean, trying to think of a way around it...
  errsum += read_config_var(v_file, keywords[36], Paramfile->pathdir);
  errsum += read_config_var(v_file, keywords[37], Paramfile->pathfile);
  //this is the file path+name of the paths file
  //XXX I THINK THIS WON"T WORK:strcat(Paramfile->pathdir,Paramfile->pathfile), also potentially look into snprintf()
  //p_file = (*(string(Paramfile->pathdir)+string(Paramfile->pathfile))).c_str(); //I think these need to be dereferenced?
  //p_file = (string(Paramfile->pathdir)+string(Paramfile->pathfile)).c_str(); //I think these need to be dereferenced?
  //p_file = (*(Paramfile->pathdir)+*(Paramfile->pathfile)).c_str(); //I think these need to be dereferenced?
  //this could work, so p_file is the parameter file. First we give it the path, and then the file that contains all the paths, paths.txt.
  strcpy(p_file,Paramfile->pathdir);
  strcat(p_file,Paramfile->pathfile);
  printf(p_file); 
  //read in the base directory
  errsum += read_config_var(p_file, keywords[38], Paramfile->basedir);
  //next up, we'll read in all the paths from the paths file

  //for sources
  errsum += read_config_var(p_file, keywords[12], Paramfile->sourcedir);
  //I think we should concatenate now... something like this
  //this puts the output into a directory named after the run? Is this what we want to do?
  //strcat(Paramfile->outputdir,Paramfile->run_name);
  //this places the base directory path before the path to the source directory. It will be repeated later for the other directories
  string(Paramfile->sourcedir).insert(0,Paramfile->basedir);
  //get the filename of the source list
  errsum += read_config_var(v_file, keywords[13], Paramfile->sourcelist);
  //because this was done after adding the basedir to the path, it should be fine?
  string(Paramfile->sourcelist).insert(0,Paramfile->sourcedir);
  //  Paramfile->sourcelist = (string(Paramfile->basedir)+string(Paramfile->sourcedir)+string(Paramfile->sourcelist)).c_str();

  //For the lenses
  errsum += read_config_var(p_file, keywords[15], Paramfile->lensdir);
  string(Paramfile->lensdir).insert(0,Paramfile->basedir);
  //Paramfile->lensdir = (string(Paramfile->basedir)+string(Paramfile->lensdir)).c_str();
  errsum += read_config_var(v_file, keywords[16], Paramfile->lenslist);
  string(Paramfile->lenslist).insert(0,Paramfile->lensdir);
  //Paramfile->lenslist = (string(Paramfile->basedir)+string(Paramfile->lensdir)+string(Paramfile->lenslist)).c_str();

  //starfields
  errsum += read_config_var(p_file, keywords[10], Paramfile->starfielddir);
  string(Paramfile->starfielddir).insert(0,Paramfile->basedir);
  //Paramfile->starfielddir = (string(Paramfile->basedir)+string(Paramfile->starfielddir)).c_str();
  errsum += read_config_var(v_file, keywords[11], Paramfile->starfieldlist);
  string(Paramfile->starfieldlist).insert(0,Paramfile->starfielddir);
  //Paramfile->starfieldlist = (string(Paramfile->basedir)+string(Paramfile->starfielddir)+string(Paramfile->starfieldlist)).c_str();

  //weather
  //XXX ADD CORRECT KEYWORD TO PATHS.TXT
  errsum += read_config_var(p_file, keywords[5], Paramfile->weatherprofiledir);
  string(Paramfile->weatherprofiledir).insert(0,Paramfile->basedir);
  //Paramfile->starfielddir = (string(Paramfile->basedir)+string(Paramfile->starfielddir)).c_str();

  //observatories
  errsum += read_config_var(v_file, keywords[0], Paramfile->obsdir);
  string(Paramfile->obsdir).insert(0,Paramfile->basedir);
  //Paramfile->obsdir = (string(Paramfile->basedir)+string(Paramfile->obsdir)).c_str();
  errsum += read_config_var(v_file, keywords[1], Paramfile->obslist);
  string(Paramfile->obslist).insert(0,Paramfile->obsdir);
  //Paramfile->obslist = (string(Paramfile->basedir)+string(Paramfile->obsdir)+string(Paramfile->obslist)).c_str();

  //planets
  

  //I need to run, but this twostep process needs to happen for all the directories in keywords[]
  //Something similar will need to be done for the actual files, as those used to be path inclusive

  //S now let's read in all the paths from the path file, for later use
  //XXX I'm 100% sure this won't work until I understand c enough
  // = (string(Paramfile->pathdir)+string(Paramfile->pathfile)).c_str()

  //WILL NEED TO CHANGE FOR NEW DIR ABSTRACTION
  //example: I think? (string(Paramfile->base_dir)+string(Paramfile->filename)).c_str()
  errsum += read_config_var(v_file, keywords[0], Paramfile->obsdir);
  errsum += read_config_var(v_file, keywords[1], Paramfile->obslist);
  if(read_config_var(v_file, keywords[2], Paramfile->setseedtoclock))
    {
      Paramfile->setseedtoclock[0]='1'; Paramfile->setseedtoclock[1]='\0';
      cerr << "Setting " << keywords[2] << " to 1 by default" << endl;
    }
  //read_config_var(v_file, keywords[3], str8); //seed -left til later
  errsum += read_config_var(v_file, keywords[4], str7); //zerotime
  //errsum += read_config_var(v_file, keywords[5], Paramfile->weatherprofiledir);
  errsum += read_config_var(v_file, keywords[6], Paramfile->run_name);
  //errsum += read_config_var(v_file, keywords[7], Paramfile->outputdir);
  if(read_config_var(v_file, keywords[8], str0)) //principle obs
    {
      str0[0]='0'; str0[1]='\0';
      cerr << "Setting " << keywords[8] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[9], str9)) //outputlc
    {
      str9[0]='0'; str9[1]='\0';
      cerr << "Setting " << keywords[9] << " to 0 by default" << endl;
    }
  //errsum += read_config_var(v_file, keywords[10], Paramfile->starfielddir);
  //errsum += read_config_var(v_file, keywords[11], Paramfile->starfieldlist);
  //errsum += read_config_var(v_file, keywords[12], Paramfile->sourcedir);
  //errsum += read_config_var(v_file, keywords[13], Paramfile->sourcelist);
  if(read_config_var(v_file, keywords[14], str1)) //source_colours
    {
      str1[0]='0'; str1[1]='\0';
      cerr << "Setting " << keywords[14] << " to 0 by default" << endl;
    }
  //errsum += read_config_var(v_file, keywords[15], Paramfile->lensdir);
  //errsum += read_config_var(v_file, keywords[16], Paramfile->lenslist);
  if(read_config_var(v_file, keywords[17], str2)) //lens_colours
    {
      str2[0]='0'; str2[1]='\0';
      cerr << "Setting " << keywords[17] << " to 0 by default" << endl;
    }
  errsum += read_config_var(v_file, keywords[18], Paramfile->planetdir); 
  errsum += read_config_var(v_file, keywords[19], Paramfile->planetroot); 
  errsum += read_config_var(v_file, keywords[20], str4); //nfilters
  errsum += read_config_var(v_file, keywords[21], str5); //Amin
  errsum += read_config_var(v_file, keywords[22], str6); //large_psf_mag
  if(read_config_var(v_file, keywords[23], str8)) //output_images
    {
      str8[0]='0'; str8[1]='\0';
      cerr << "Setting " << keywords[23] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[24], str10)) //prettypic
    {
      str10[0]='0'; str10[1]='\0';
      cerr << "Setting " << keywords[24] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[25], str11)) //prettypicDim
    {
      strcpy(str11,"256");
      cerr << "Setting " << keywords[25] << " to 256 by default" << endl;
    }
  errsum += read_config_var(v_file, keywords[26], str12); //min_chi2
  if(read_config_var(v_file, keywords[27], str13)) //outputOnErr
    {
      str13[0]='0'; str13[1]='\0';
      cerr << "Setting " << keywords[27] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[28], str14)) //outputOnDet
    {
      str14[0]='0'; str14[1]='\0';
      cerr << "Setting " << keywords[28] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[29], str15)) //outputOnAll
    {
      str15[0]='0'; str15[1]='\0';
      cerr << "Setting " << keywords[29] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[30], str17)) //pllxMultiplyer
    {
      str17[0]='1'; str17[1]='\0';
      cerr << "Setting " << keywords[30] << " to 1 by default" << endl;
    }
  if(read_config_var(v_file, keywords[31], str16)) //lenslight
    {
      str16[0]='1'; str16[1]='\0';
      cerr << "Setting " << keywords[31] << " to 1 by default" << endl;
    }
  if(read_config_var(v_file, keywords[32], str18)) //repeat sequence
    {
      str18[0]='0'; str18[1]='\0';
      cerr << "Setting " << keywords[32] << " to 0 by default" << endl;
    }
  if(read_config_var(v_file, keywords[33], Paramfile->obsgroupstr))
    {
      strcpy(Paramfile->obsgroupstr,"(ALL)");
    }
  if(read_config_var(v_file, keywords[34], str19)) //NUM_SIM_DAYS
    {
      strcpy(str19,"2010");
    }
  if(read_config_var(v_file, keywords[35], str20)) //U0MAX
    {
      strcpy(str20,"3");
    }


  /*char *keywords[27] = {"OBSERVATORY_DIR", "OBSERVATORY_LIST", 
			"SET_RANDOM_SEED_TO_CLOCK", "RANDOM_SEED", 
			"SIMULATION_ZERO_TIME", "WEATHER_PROFILE_DIR", 
			"RUN_NAME", "OUTPUT_DIR", "PRINCIPLE_OBSERVATORY",
			"OUTPUT_LC", "STARFIELD_DIR", "STARFIELD_LIST",
			"SOURCE_DIR", "SOURCE_LIST", "SOURCE_COLOURS", 
			"LENS_DIR", "LENS_LIST", "LENS_COLOURS", "PLANET_DIR",
			"PLANET_ROOT", "NFILTERS", "AMIN", "LARGEPSFMAG",
			"OUTPUT_IMAGES", "PRETTY_PICS", 
			"PRETTY_PICS_DIMENSIONS", "MIN_CHISQUARED",
			"OUTPUT_ONERR", "OUTPUT_ONDET", "OUTPUT_ONALL"
			"PARALLAX", "LENS_LIGHT","REPEAT_SEQUENCE",
			"OBS_GROUPS","NUM_SIM_DAYS"};*/

  //do any additional processing

  strcat(Paramfile->outputdir,Paramfile->run_name);
  strcat(Paramfile->outputdir,"/");
  sscanf(Paramfile->setseedtoclock,"%d",&Paramfile->setseedtoclockBIT);
  sscanf(str7, "%lf", &Paramfile->simulation_zerotime);
  Paramfile->outputLightcurve=atof(str9);
  Paramfile->principle_observatory=atoi(str0);
  Paramfile->sourcecolours=atoi(str1);
  Paramfile->lenscolours=atoi(str2);
  Paramfile->Nfilters = atoi(str4);
  Paramfile->Amin = atof(str5);
  Paramfile->large_psf_mag = atof(str6);
  Paramfile->outputImages = atoi(str8);
  Paramfile->prettypic = atoi(str10);
  Paramfile->minChiSquared = atof(str12);
  Paramfile->outputOnErr = atoi(str13);
  Paramfile->outputOnDet = atoi(str14);
  Paramfile->outputOnAll = atoi(str15);
  Paramfile->pllxMultiplyer = atoi(str17);
  Paramfile->lenslight = atoi(str16);
  Paramfile->identicalSequence = atoi(str18);
  Paramfile->NUM_SIM_DAYS = atoi(str19);
  Paramfile->u0max = atof(str20);

  string string11 = string(str11);
  size_t pos;
  pos = string11.find_first_of(",xX:");
  if(pos==string::npos)
    {
      Paramfile->prettypicDimX = atoi(str11);
      Paramfile->prettypicDimY = Paramfile->prettypicDimX;
    }
  else
    {
      Paramfile->prettypicDimX = atoi(string11.substr(0,pos).c_str());
      Paramfile->prettypicDimY = atoi(string11.substr(pos+1).c_str());
    }

  if(Paramfile->setseedtoclockBIT==0)
    {
      read_config_var(v_file, keywords[3] , str3);
      sscanf(str3,"%ld",&Paramfile->Seed);
    }

  sprintf(str,"Input file %s",v_file);
  fmtline(str,WIDTH,"PARSED");  

  if(errsum<0)
    {
      cerr << "There were essential parameters missing. Exiting." << endl;
      exit(1);
    }

}
