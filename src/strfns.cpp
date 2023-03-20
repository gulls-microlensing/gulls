#include "strfns.h"

#include<iostream>

using namespace std;


char *pad(char *s, int size)
{
   int len = strlen(s);
   if (s != NULL && len < size)
      sprintf(s+len,"%*c",size-len,' ');
   return s;
}

char *fmtline(char *instr, int size, const char *msg){

 pad(instr,size);
 printf(strcat(strcat(instr,msg),"\n"));
 return(0);

}


int read_config_var( char *values_file, const char *keyword , char value[] )
{
  static char str[1000];
  int len;
  FILE * _file=NULL; 

 

  if( keyword == NULL ) 
    {
      cerr << "read_config_var: Error: No keyword passed" << endl;
      return (KEY_ERR);
    }

  len = strlen(keyword);

  if( len > 77)
    {
      cerr << "read_config_var: Error: Keyword "  << keyword << " too long" << endl;
      return(KEY_ERR);
    }


  if( values_file )
    {   
      _file = fopen(values_file, "r");
      if (_file == NULL)
	{ 
	  cerr << "read_config_var: Error: Could not open parameter file" << endl;
	  return(FILE_ERR);
	}
    }


  if( fseek(_file, 0, SEEK_SET) )
    {
      fclose(_file);
      cerr << "read_config_var: Error: Unknown file seek error" << endl;
      return(FILE_ERR);
    }

  for(;;)
    {
      fgets(str, 1000, _file);
      if( ferror(_file) )
	{
	  cerr << "read_config_var: Error: Error reading file" << endl;
	  return(FILE_ERR);
	}
      if(feof(_file))
	{
	  cerr << "read_config_var: Error: Keyword " << keyword << " not found" << endl;
	  return(FILE_ERR);
	}
      len = strlen(str);
      
      if( strncmp(keyword, str, strlen(keyword)) == 0 )
	{ 
	  if (str[len - 1] == '\n') str[--len] = 0;
	  sprintf(value, "%s", &str[strlen(keyword)+1] ); break;
	}
    } 
  
  fclose(_file);
  
  return 0;

}

int strip_whitespace(char* input)
{
  char str[1000];
  int i=0;
  int j=0;
  
  strcpy(str,input);

  /*strip out whitespace \t <space>*/
  while(str[i]!=0)
    {
      if(str[i]!=' ' && str[i]!='\t')
	{
	  input[j]=str[i];
	  j++;
	} 
      i++;
    }

  input[j]=0;

  return 0;

}
