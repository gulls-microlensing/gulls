#ifndef COLUMNCODES

//For the Galaxy model input:
static const int NDATAFIELDS = 26;

//The meanings of each of the input and output columns
static const int MUL = 0;
static const int MUB = 1;
static const int VR = 2;
static const int UU = 3;
static const int VV = 4;
static const int WW = 5;
static const int MV = 6;
static const int CL = 7;
static const int TYP = 8;
static const int TEFF = 9;
static const int LOGG = 10;
static const int AGE = 11;
static const int MASS = 12;
static const int MBOL = 13;
static const int RADIUS = 14;
static const int FEH = 15;
static const int LL = 16;
static const int BB = 17;
static const int RA = 18;
static const int DEC = 19;
static const int DIST = 20;
static const int XX = 21;
static const int YY = 22;
static const int ZZ = 23;
static const int AV = 24;
static const int HEFE = 25;

static const int sOutputCols=7;
static const int lOutputCols=11;

const int sOutputColumns[sOutputCols] = {DIST, RADIUS, MUL, MUB, AGE, CL, TYP};
const int lOutputColumns[lOutputCols] = {DIST, MASS, MUL, MUB, AGE, CL, TYP, MBOL, TEFF, LOGG, RADIUS};

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
