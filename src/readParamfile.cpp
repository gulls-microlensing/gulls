#include "readParamfile.h"


void readParamfile(char *v_file, struct filekeywords *Paramfile){
  //for reading and later combing keywords
  char temp[1000];

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
  char str21[20];
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
			      "BASE_DIR","USE_FIELDS","ERROR_SCALING"};          //36, 37, 38



  int errsum=0; //use this to check that critical keywords are included in the parameter file
                //test read_config_var return to set defaults if the keyword is not

  //getenv here for the path to get base_path, tack on to the beginning of everything, should not be bad
  //this way, we can ditch paths.txt or whatever, will make everything simpler. 


  //start off by getting the location of paths.txt, where are the directory/path info is kept
  //from here on, directories are read from p_file, while v_file is more specific
  //errsum += read_config_var(v_file, keywords[36], Paramfile->pathdir);
  //errsum += read_config_var(v_file, keywords[37], Paramfile->pathfile);
  //Paramfile is a pointer to a structure
  //printf(Paramfile->pathdir);
  //printf("%s",Paramfile->pathdir);
  //strcpy(p_file,Paramfile->pathdir);
  //strcat(p_file,Paramfile->pathfile);
  //read the parameters from the paths file paths.txt
  //errsum += read_config_var(p_file, keywords[38], Paramfile->basedir);
  //errsum += read_config_var(p_file, keywords[38], temp);
  //strcat(Paramfile->basedir,temp);
  //printf(Paramfile->basedir);
  //printf("\n");

  //get the base directory from the env variable set in .gulls
  

  //if(!std::getenv("GULLS_BASE_DIR")){
  // printf("HERE\n");
  //   }


  Paramfile->basedir = std::getenv("GULLS_BASE_DIR");

  //okay, so this is the best way i have come up with to combine info from
  //the p_file and v_file. it's not very elegant, but seems to be working okay.
  //I should write a separate function, but for as few a times as we need this work 
  //around, I am just going to write it in.


  //for observatories
  //[0] is obsdir
  //errsum += read_config_var(v_file, keywords[0], Paramfile->obsdir);
  errsum += read_config_var(v_file, keywords[0], temp);
  //need to first cat on the base directory
  strcat(Paramfile->obsdir,Paramfile->basedir);
  //then we cat on the actual directory
  strcat(Paramfile->obsdir,temp);


  //[1] is obslist
  //errsum += read_config_var(v_file, keywords[1], Paramfile->obslist);
  errsum += read_config_var(v_file, keywords[1], temp);
  strcat(Paramfile->obslist,Paramfile->obsdir);
  strcat(Paramfile->obslist,temp);

  if(read_config_var(v_file, keywords[2], Paramfile->setseedtoclock))
    {
      Paramfile->setseedtoclock[0]='1'; Paramfile->setseedtoclock[1]='\0';
      cerr << "Setting " << keywords[2] << " to 1 by default" << endl;
    }

  //read_config_var(v_file, keywords[3], str8); //seed -left til later
  errsum += read_config_var(v_file, keywords[4], str7); //zerotime
  //[5] is weatherdir
  //errsum += read_config_var(v_file, keywords[5], Paramfile->weatherprofiledir);
  errsum += read_config_var(v_file, keywords[5], temp);
  strcat(Paramfile->weatherprofiledir,Paramfile->basedir);
  strcat(Paramfile->weatherprofiledir,temp);

  errsum += read_config_var(v_file, keywords[6], Paramfile->run_name);
  //left output dir separate from restructuring for now
  errsum += read_config_var(v_file, keywords[7], Paramfile->outputdir);
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

  //[10] is starfielddir
  //errsum += read_config_var(v_file, keywords[10], Paramfile->starfielddir);
  errsum += read_config_var(v_file, keywords[10], temp);
  strcat(Paramfile->starfielddir,Paramfile->basedir);
  strcat(Paramfile->starfielddir,temp);

  //[11] is starfieldlist
  //errsum += read_config_var(v_file, keywords[11], Paramfile->starfieldlist);
  errsum += read_config_var(v_file, keywords[11], temp);
  strcat(Paramfile->starfieldlist,Paramfile->starfielddir);
  strcat(Paramfile->starfieldlist,temp);

  //[12] is sourcedir
  //errsum += read_config_var(v_file, keywords[12], Paramfile->sourcedir);
  errsum += read_config_var(v_file, keywords[12], temp);
  strcat(Paramfile->sourcedir,Paramfile->basedir);
  strcat(Paramfile->sourcedir,temp);

  //[13] is sourcelist
  //errsum += read_config_var(v_file, keywords[13], Paramfile->sourcelist);
  errsum += read_config_var(v_file, keywords[13], temp);
  strcat(Paramfile->sourcelist,Paramfile->sourcedir);
  strcat(Paramfile->sourcelist,temp);

  if(read_config_var(v_file, keywords[14], str1)) //source_colours
    {
      str1[0]='0'; str1[1]='\0';
      cerr << "Setting " << keywords[14] << " to 0 by default" << endl;
    }


  //[15] is lensdir
  //errsum += read_config_var(v_file, keywords[15], Paramfile->lensdir);
  errsum += read_config_var(v_file, keywords[15], temp);
  strcat(Paramfile->lensdir,Paramfile->basedir);
  strcat(Paramfile->lensdir,temp);

  //[16] is lenslist
  //errsum += read_config_var(v_file, keywords[16], Paramfile->lenslist);
  errsum += read_config_var(v_file, keywords[16], temp);
  strcat(Paramfile->lenslist,Paramfile->lensdir);
  strcat(Paramfile->lenslist,temp);

  if(read_config_var(v_file, keywords[17], str2)) //lens_colours
    {
      str2[0]='0'; str2[1]='\0';
      cerr << "Setting " << keywords[17] << " to 0 by default" << endl;
    }


  //[18] is planetdir
  //errsum += read_config_var(v_file, keywords[18], Paramfile->planetdir);
  errsum += read_config_var(v_file, keywords[18], temp);
  strcat(Paramfile->planetdir,Paramfile->basedir);
  strcat(Paramfile->planetdir,temp);

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

  if(read_config_var(v_file, keywords[38], str21)) //ERROR_SCALING
    {
      strcpy(str21,"0");
      cerr << "Setting " << keywords[38] << " to 0 by default" << endl;

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
  Paramfile->error_scaling = atoi(str21);

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
