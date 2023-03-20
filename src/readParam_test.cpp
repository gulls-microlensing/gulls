#include <stdio.h>
#include "strfns.h"
#include "structures.h"
#include "readParamfile.h"

int main(){
  //just going to test readParamfile.c
  printf("Hello makefiles!\n");
  //std::string input_filename; 
  char *input_filename;
  const char *env_test = getenv("GULLS_BASE_DIR");
  printf("env test\n");
  printf(env_test);
  printf("env test comlpete\n");
  printf("\n");
  input_filename = "/home/stingray/johnson.7080/gulls/gulls_sj/params/parameterFiles/kmtwfffp.prm";    /* FILE HANDLING */
  //printf(input_filename);
  //printf("\n");
  struct filekeywords Paramfile;
  readParamfile(input_filename, &Paramfile);
  printf("POST READ PARAM: \n");
  printf(Paramfile.pathdir);
  printf("\n");
  printf(Paramfile.pathfile);
  printf("\n");
  return 0;
}
