#ifndef COLUMNCODES

//For the Galaxy model input:
//static const int NDATAFIELDS = 26;



//static const int sOutputCols=7;
//static const int lOutputCols=11;

//const int sOutputColumns[sOutputCols] = {DIST, RADIUS, MUL, MUB, AGE, CL, TYP};
//const int lOutputColumns[lOutputCols] = {DIST, MASS, MUL, MUB, AGE, CL, TYP, MBOL, TEFF, LOGG, RADIUS};

//Planet data input
static const int NPLANETINPUT = 4;
static const int NPLANETDERIV = 3;

//input
static const int PMASS = 0;
static const int AA = 1; //semimajor axis
static const int INC = 2;
static const int PHASE = 3;
//derived
static const int QQ = NPLANETINPUT + 0; //mass ratio
static const int SS = NPLANETINPUT + 1; //separation
static const int TT = NPLANETINPUT + 2; //Period


#define COLUMNCODES
#endif
