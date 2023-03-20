#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef NOERROR
	#define NOERROR 0
#endif
#ifndef FILE_ERR
	#define  FILE_ERR -1
#endif
#ifndef KEY_ERR
	#define KEY_ERR -2
#endif


char *pad(char *s, int size);
char *fmtline(char *str, int size, const char *msg);
int read_config_var( char *, const char * , char [] );
int strip_whitespace(char* input);
