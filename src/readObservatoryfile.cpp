#include "readObservatoryfile.h"
#include "constdefs.h"
#include<iostream>

#define DEBUGVAR 0

using namespace std;

void specerror(int errval,const char keyword[]){

  switch(errval){
 case 1:
   printf("(readObservatoryfile): Keyword error: keyword %s not found. Exiting\n",keyword);
  break;
 case 2:
   printf("(readObservatoryfile): Keyword error: %s conversion failure. Exiting\n",keyword);
  break;
  }
}

void readObservatoryfile(string v_file, struct obsfilekeywords World[],int idx){
  const char *keywords[26] = {"NAME","LATITUDE","LONGITUDE","ALTITUDE","READ_OHEAD","NFIELDS","WEATHER_PROFILE","FIELDCENTRES","NPIX_X","NPIX_Y","PIXELSIZE","PRIMARY","BLOCKAGE","SPACE","FILTER","OBSERVATION_SEQUENCE","DETECTOR","THROUGHPUT","REFERENCE_TEXP","REFERENCE_NSTACK","ORBIT","PHOTOMETRY","EXTCOEFF","SKY_BACKGROUND","MOON_AVOID","ALT_LIMIT"};

  int nkey=26; /* Number of keywords defined in array "keywords" */
  char String[26][1000];      /*!< Matrix of string variables */

  int jdx;
	  
  /* Read the Keyword values in as strings */
  for(jdx=0;jdx<nkey;jdx++)
    {
      char v_file_temp[1000];
      if(read_config_var(v_file.c_str(), keywords[jdx] , String[jdx])!=0)
	{
	  cerr << "Error reading observatory file (" << v_file << ")" << endl;
	  specerror(1,keywords[jdx]);
	  switch(jdx)
	    {
	    case 22:
	    case 24:
	    case 25:
	      break;
	    default:
	      exit(1);
	    }

	}
    }

  /*  Have to convert some to doubles in the observatory structure */

  World[idx].name = string(String[0]); //"NAME"
  if(sscanf(String[1],"%lf",&World[idx].latitude) ==0) {specerror(2,keywords[1]); exit(1);}; //"LATITUDE"
  if(sscanf(String[2],"%lf",&World[idx].longitude) ==0) {specerror(2,keywords[2]); exit(1);}; //"LONGITUDE"
  if(sscanf(String[3],"%lf",&World[idx].altitude) ==0) {specerror(2,keywords[3]); exit(1);}; //"ALTITUDE"

  if(sscanf(String[4],"%lf",&World[idx].readohead) ==0) {specerror(2,keywords[4]); exit(1);}; //"READ_OHEAD"

  if(sscanf(String[5],"%d",&World[idx].nfields) ==0) {specerror(2,keywords[5]); exit(1);}; //"NFIELDS"    
  World[idx].weatherProfile = string(String[6]); //"WEATHER_PROFILE"

  //"NPIX_X","NPIX_Y","PIXELSIZE","PRIMARY","BLOCKAGE","SPACE","FILTER","OBSERVATION_SEQUENCE","DETECTOR"

  if(sscanf(String[8],"%d",&World[idx].npixx) ==0) {specerror(2,keywords[8]); exit(1);}; //"NPIX_X"
  if(sscanf(String[9],"%d",&World[idx].npixy) ==0) {specerror(2,keywords[9]); exit(1);}; //"NPIX_Y"
  if(sscanf(String[10],"%lf",&World[idx].pixelsize) ==0) {specerror(2,keywords[10]); exit(1);}; //"PIXELSIZE"

  if(sscanf(String[11],"%lf",&World[idx].primary) ==0) {specerror(2,keywords[11]); exit(1);}; //"PRIMARY"
  if(sscanf(String[12],"%lf",&World[idx].blockage) ==0) {specerror(2,keywords[12]); exit(1);}; //"BLOCKAGE"

  if(sscanf(String[13],"%d",&World[idx].space) ==0) {specerror(2,keywords[13]); exit(1);}; //"SPACE"

  if(sscanf(String[14],"%d",&World[idx].filter)==0) {specerror(2,keywords[14]); exit(1);}; //"FILTER"

  World[idx].observationSequence = string(String[15]); 
  //"OBSERVATION_SEQUENCE"
  World[idx].detector = string(String[16]); //"DETECTOR"
  World[idx].throughput = string(String[17]); //"THROUGHPUT"

  if(sscanf(String[18],"%lf",&World[idx].reftexp) ==0) {specerror(2,keywords[18]); exit(1);}; //"REFERENCE_TEXP"
  if(sscanf(String[19],"%lf",&World[idx].refnstack) ==0) {specerror(2,keywords[19]); exit(1);}; //"REFERENCE_NSTACK"

  World[idx].orbitcode = string(String[20]); //ORBIT

  if(sscanf(String[21],"%d",&World[idx].photcode)==0) {specerror(2,keywords[14]); exit(1);}; //"PHOTOMETRY"
  if(sscanf(String[22],"%lf",&World[idx].extcoeff) ==0) {specerror(2,keywords[22]); World[idx].extcoeff=0.0;}; //"EXTCOEFF"
  if(sscanf(String[23],"%lf",&World[idx].skybackground) ==0) {specerror(2,keywords[23]); World[idx].skybackground=99.0;}; //"SKY_BACKGROUND"
  if(sscanf(String[24],"%lf",&World[idx].moonavoid) ==0) {specerror(2,keywords[24]); World[idx].moonavoid=30.0;}; //"MOON_AVOID"
    if(sscanf(String[25],"%lf",&World[idx].altlimit) ==0) {specerror(2,keywords[25]); World[idx].altlimit=30.0;}; //"MOON_AVOID"

}




