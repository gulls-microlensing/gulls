#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "structures.h"

using namespace std;

int readSLList(int choosefield, bool src, struct filekeywords *Paramfile, struct slcat *sl);
int getValidFields(struct filekeywords* Paramfile, struct slcat *Sources, struct slcat *Lenses, vector<double>* sfielddata);
