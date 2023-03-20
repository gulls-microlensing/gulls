#include "readSynthGal.h"
#include<math.h>

#define DEBUGVAR 0

int readSynthGal(struct filekeywords *Paramfile, struct galaxy *Galaxy, char* instance){
  int file_exists(const char * filename);
  FILE *istream,*istreamlock;
  char galfracfilename[1000];
  //char galfracfilename1[1000];

  char galfracfilenamelock[1000];

  int i;
  //char suffixstr[10];
  //int maxSubGalaxies;
  //double tmp;

  strcpy(galfracfilename,Paramfile->synthgaldir);
  strcat(galfracfilename,Paramfile->synthgalroot);
  /* strcat(galfracfilename,"."); */
  strcat(galfracfilename,instance);

  if(DEBUGVAR) printf("Galaxy filename: %s\n",galfracfilename);

  /*  TEMPORARILY DISABLING FILE LOCKING AS IS NOT WORKING PROPERLY */
  /*  GALAXY SUFFIX PASSED DIRECTLY AS INSTANCE                     */
  /*
  tmp=pow(10,Paramfile->suffixLength); 
  maxSubGalaxies = (int)tmp;

  for(suffix=0;suffix<maxSubGalaxies;suffix++){
    sprintf(suffixstr,Paramfile->suffixForm,suffix);
  strcpy(galfracfilename1,galfracfilename);

  strcat(galfracfilename1,suffixstr);                               */
  /*  printf("reading galaxy suffix %d format %s name %s\n",suffix,Paramfile->suffixForm,galfracfilename1); */
  /*  sprintf(galfracfilenamelock,"%s.lock",galfracfilename1);      */
  sprintf(galfracfilenamelock,"%s.lock",galfracfilename);


  /*  if ( (istream = fopen ( galfracfilename1, "r" ) ) != NULL ){*/
  if ( (istream = fopen ( galfracfilename, "r" ) ) != NULL )
    {
      /* printf("%s Exists\n", galfracfilename1); */
      if ( (istreamlock = fopen ( galfracfilenamelock, "r" ) ) == NULL )
	{
	  /*  printf("%s Does not exist\n", galfracfilenamelock); */
	  printf("Reading file: %s\n",galfracfilename);
	  istreamlock=fopen( galfracfilenamelock, "w" );
	  fclose(istreamlock);
	  
	  i=0;
	  while(fscanf(istream,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %d\n", &Galaxy->l[i], &Galaxy->b[i], &Galaxy->dS[i], &Galaxy->dL[i], &Galaxy->m_S[i].m[0],  &Galaxy->m_S[i].m[1], &Galaxy->m_S[i].m[2], &Galaxy->m_S[i].m[3], &Galaxy->m_S[i].m[4], &Galaxy->m_L[i].m[0], &Galaxy->m_L[i].m[1], &Galaxy->m_L[i].m[2], &Galaxy->m_L[i].m[3], &Galaxy->m_L[i].m[4], &Galaxy->Mbol_S[i], &Galaxy->Rs[i], &Galaxy->mL[i], &Galaxy->vperp[i], &Galaxy->tE[i], &Galaxy->umax[i], &Galaxy->w[i], &Galaxy->id[i], &Galaxy->mp[i], &Galaxy->a[i], &Galaxy->inc[i], &Galaxy->phase[i], &Galaxy->eA[i], &Galaxy->fieldno[i]) !=EOF)
	    {
	      if(i>MAXGALVALS) 
		{
		  printf("(readSynthGal) Array overload (MAXGALVALS=%d)\n", MAXGALVALS); 
		  exit(1);
		}
	      i++;
	    }
	  fclose(istream);

	  Galaxy->nsynth = i;
	  /*  printf("Number of galaxy entries: %d\n",Galaxy->nsynth); */
	  printf("Galaxy read in\n"); fflush(stdout);
	  return(1); /* i.e. do not continue loop over suffix values */
	}  
       else printf("%s Exists\n",galfracfilenamelock); 
    }
    else printf("%s file does not exist\n", galfracfilename);

    /*  }*/

  return(0);
}
