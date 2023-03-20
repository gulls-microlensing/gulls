/*! \file 
\brief A collection of astronomy-related functions

This is a repository of all the main astronomy functions
 */

#include "astroFns.h"
#include <gsl/gsl_blas.h>
#include "integerPowers.h"

/*! Compute Local Sidereal Time */
void lst(double JD, double EastLong, double *LST){

double TJD,DayFrac,T,GMST0UT;

/* -------------------------------------------------------------------- */
/*  lst function         Local Sidereal Time, (mean or apparent), */
/*                     for vector of JD's and a given East Longitude. */
/*  input  : - Vector of JD, in UT1 time scale. */
/*           - East Longitude in radians. */
/*           - Sidereal Time Type, */
/*             'm' - Mean (default). */
/*             'a' - apparent. */
/*  output : - vector of LST in fraction of day. */
/* -------------------------------------------------------------------- */

/*  convert JD to integer day + fraction of day */
TJD = floor(JD - 0.5) + 0.5;
DayFrac = JD - TJD;

T = (TJD - 2451545.0)/36525.0;

GMST0UT = 24110.54841 + 8640184.812866*T + 0.093104*T*T - 6.2e-6*T*T*T;

 /* convert to fraction of day in range [0 1) */
GMST0UT = GMST0UT/86400.0;

GMST0UT = GMST0UT - floor(GMST0UT);
*LST = GMST0UT + 1.0027379093*DayFrac + EastLong/(2*PI);
*LST = *LST - floor(*LST);

}



void eq2horiz(double incoord1, double incoord2, double JD, double Lat, double Long, char Direction, double *outcoord1, double *outcoord2){

/* -------------------------------------------------------------------- */
/*  eq2horiz function      Horizontal coordinates conversion */
/*                        Converting from equatorial coordinates */
/*                        to Horizontal coordinates and visa */
/*                        versa. */
/*  input  : -  coordinates: (RA & Dec) | (Az & Alt) in radians */
/*           - JD + UT fraction of day, */
/*           - Geodetic Coordinates, east long & north lat in radians */
/*           - Direction, */
/*             'h' - from equatorial to horizontal (default). */
/*             'e' - from horizontal to equatorial. */
/*  output : - converted coordinates. */
/* Adapted from horiz_coo.m */
/* -------------------------------------------------------------------- */

  double LST,HA,RA,Dec,Alt,Az,SinAlt,CosAlt,SinAz,CosAz;
  double SinDec,CosDec, SinHA, CosHA;

/*  calculating Local Mean Sidereal Time */
  lst(JD,Long,&LST);


switch(Direction){
 case 'h':

    /* convert equatorial to horizontal */

    /* calculate the Hour Angle */
   HA  = LST*2*PI - incoord1;
   Dec = incoord2;


   SinAlt = sin(Dec)*sin(Lat) + cos(Dec)*cos(HA)*cos(Lat);
   CosAlt = sqrt(1-SinAlt*SinAlt);

   SinAz  = (-cos(Dec)*sin(HA))/CosAlt;
   CosAz  = (sin(Dec)*cos(Lat) - cos(Dec)*cos(HA)*sin(Lat))/CosAlt;

   Az     = atan2(SinAz, CosAz);
   Alt    = asin(SinAlt);

   while(Az<0){
     Az+=2*PI;
   }

   *outcoord1 = Az;
   *outcoord2 = Alt;

   break;

 case 'e':

   Az     = incoord1;
   Alt    = incoord2;

   SinDec = sin(Alt)*sin(Lat) + cos(Alt)*cos(Az)*cos(Lat);
   CosDec = sqrt(1 - SinDec*SinDec);

   SinHA  = (-cos(Alt)*sin(Az))/CosDec;
   CosHA  = (sin(Alt)*cos(Lat) - cos(Alt)*cos(Az)*sin(Lat))/CosDec;
   HA     = atan2(SinHA, CosHA);

   RA     = LST*2*PI - HA;
   Dec    = asin(SinDec);

    /* converting to range [0,1) */
   RA     = 2*PI*(RA/(2*PI) - floor(RA/(2*PI)));

   *outcoord1 = RA;
    *outcoord2= Dec;
    break;

 default:
   printf("Illigal type of conversion\n"); exit(1);
   break;
  }

}

void eq2gal(double in_coord1, double in_coord2, char dir, double *out_coord1, double *out_coord2)
{

  /*  Input: (RA & DEC) || (GAL_LONG & GAL_LAT) in radians */
  /*     Output: converted values in radians */
  /*     dir = 'e' for Equatorial to Galactic; 'g' for Galactic to Equatorial*/
  /*  NOTE: eq2gal.m is WRONG - it has these the wrong way around */
 

  double a[] = {
    -0.0548755604,+0.4941094279,-0.8676661490,
    -0.8734370902,-0.4448296300,-0.1980763734,
    -0.4838350155,+0.7469822445,+0.4559837762};

  double b[3];
  double ans[] = {0.0, 0.0, 0.0};
  
  gsl_matrix_view Mat = gsl_matrix_view_array(a, 3, 3);

  b[0] = cos(in_coord1)*cos(in_coord2);
  b[1] = sin(in_coord1)*cos(in_coord2);
  b[2] = sin(in_coord2);

  gsl_vector_view Pos = gsl_vector_view_array(b, 3);
  
  gsl_vector_view RPos = gsl_vector_view_array(ans, 3);


  /* printf("%f %f %f\n",b[0],b[1],b[2]); */

  switch(dir){
  case 'g':
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Mat.matrix, &Pos.vector, 1.0, &RPos.vector);
    break;
  case 'e':
    gsl_blas_dgemv(CblasTrans, 1.0, &Mat.matrix, &Pos.vector, 1.0, &RPos.vector);
    break;
  default:
    printf("unknown type of conversion\n");
    break;
  }

  /*  printf ("[%g %g %g]\n", ans[0], ans[1], ans[2]); */

     
  /* % convert to long & lat */
  /* RPos = RPos'; */
  
  *out_coord1 = atan2(ans[1],ans[0]);
  *out_coord2 = atan(ans[2]/sqrt(ans[0]*ans[0] + ans[1]*ans[1]));
 
}


void suncoo(double JD, char EquinoxType, double *RA, double *Dec, double *R,double *SL, double *EquationTime){

/* -------------------------------------------------------------------- */
/*  suncoo function   Calculate Sun Equatorial coordinates */
/*                  low accuracy formulae. */
/*                  Accuracy : 0.01 deg. in long. */
/*  input  : - Vector of JDs. */
/*           - Equinox for output coordinates: */
/*             'g' - True Geometric referred to the mean equinox of date. */
/*             'a' - Apparent, referred to the true equinox of date. (default). */
/*             'j' - J2000, referred to the J2000.0 equinox. */
/*  output : - vector of RA, in radians. */
/*           - vector of Dec. in radians. */
/*           - Vector of radius vectors, in AU. */
/*           - Solar longitude in the same ref. frame as RA/Dec. (radians). */
/*           - Equation of Time [Minuts of time] */
/* -------------------------------------------------------------------- */

  double T,L0,M,e,C,Ni,Om,Obl;

T   = (JD - 2451545.0)/36525;
 L0  = (280.46645 + 36000.76983*T + 0.0003032*T*T)/(TO_DEG);
 M   = (357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T)/(TO_DEG);
e   = 0.016708617 - 0.000042037*T - 0.0000001236*T*T;


C   = (1.914600 - 0.004817*T - 0.000014*T*T)*sin(M) + (0.019993 - 0.000101*T)*sin(2*M) + 0.000290*sin(3*M);
C   = C/(TO_DEG);
 /* Sun longitude */
*SL  = L0 + C;

 /* the sun true Anomaly: */
Ni  = M + C;

/*  the sun radius vector */
*R   = 1.000001018*(1-e*e)/(1+e*cos(Ni));

switch(EquinoxType){
 case 'a':
   Om = (125.04 - 1934.136*T)/(TO_DEG);
   *SL = *SL - (0.00569 - 0.00478*sin(Om))/(TO_DEG);
   break;
 case 'j':
   *SL = *SL - (0.01397*T*100)/(TO_DEG);
   break;
 case 'g':
    /* Already geometric longitude */
   break;
 default:
   printf("Illegal equinox type\n"); exit(1);
   break;
  }


 Obl = obliquity(JD,'L');

 *SL     = (*SL/(2*(PI)) - floor(*SL/(2*(PI))))*2*PI;
*RA     = atan2(cos(Obl)*sin(*SL),cos(*SL));
*Dec    = asin(sin(Obl)*sin(*SL));


*EquationTime = *RA - *SL;

 while(*EquationTime < -PI) { *EquationTime += 2*PI;}
 while(*EquationTime > PI) { *EquationTime -= 2*PI;}

*EquationTime *= 360;

}

double obliquity(double JulianDay, char Type){

/* --------------------------------------------------------- */
/*  obliquity function      calculating obliquity of */
/*                          ecliptic (with respect to the */
/*                          mean equator od date) */
/*                          for a given julian day. */
/*  Input  : - vector of Julian Days. */
/*           - Caqlculation type: */
/*             'L' - IAU 1976, good from 1000-3000 AD, */
/*                   default. */
/*             'H' - Laskar expression, more accurate. */
/*  Output : - obliquity of ecliptic of date in radians. */
/* --------------------------------------------------------- */

  double T,Obl,U;

switch(Type){
 case 'L':
    T   = (JulianDay - 2451545.0)/36525.0;
    //Obl = 23.439291 - 0.0130042*T - 0.00000016*T*T + 0.000000504*T*T*T;
    Obl = 23.439291 + T*(-0.0130042 + T*(-0.00000016 + 0.000000504*T));
    Obl = Obl/(TO_DEG);
    break;

 case 'H':
    T   = (JulianDay - 2451545.0)/36525.0;
    U   = T/100;
    //Obl = 23.44484666666667 - (4680.93*U - 1.55*U*U + 1999.25*pow(U,3) - 51.38*pow(U,4) - 249.67*pow(U,5) - 39.05*pow(U,6) + 7.12*pow(U,7) + 27.87*pow(U,8) + 5.79*pow(U,9)  + 2.45*pow(U,10))/3600;
    Obl = 23.44484666666667 - U*(4680.93 + U*(-1.55 + U*(1999.25 + U*(-51.38 + U*(-249.67 + U*(-39.05 + U*(7.12 + U*(27.87 + U*(5.79  + 2.45*U)))))))))/3600;
    Obl = Obl/(TO_DEG);
    break;
 default:

    printf("Unknown calculation type in obliquity\n");
    exit(1);
    break;

 }  

  return(Obl);
}



void jd2date(double JD, int *day, int *month, int *year, double *frac)
{

/* %-------------------------------------------------------------------- */
/* % jd2date function       convert Julian Day to date. */
/* % input  : - vector of Julian Days. */
/* % output : - matrix in which the first column is the Day in the month */
/* %            the second column is the month, the third is the year */
/* %            and the fourth column is fraction of day (U.T.) */
/* %    By  Eran O. Ofek           September 1999 */
/* %-------------------------------------------------------------------- */

  double Z,F,A,Day,Month,Year;
  double Alpha,B,C,D,E;

  if (JD<0)
    {
      printf("The method i valid only for poitive JDs\n");
      exit(1);
    }

  Z = floor(JD+0.5);
  F = (JD+0.5) - floor(JD+0.5);

  if(Z<2299161)  A = Z;
  else
    {
      Alpha = trunc((Z - 1867216.25)/36524.25);
      A = Z + 1 + Alpha - trunc(Alpha/4);
    }


  B = A + 1524;
  C = trunc((B - 122.1)/365.25);
  D = trunc(365.25*C);
  E = trunc((B-D)/30.6001);

  Day   = B - D - trunc(30.6001*E) + F;

  if(E<14) Month = E - 1;
  else Month = E - 13;

  if(Month>2) Year = C - 4716;
  else Year = C - 4715; //if month==1 or ==2


  *day = int(floor(Day));
  *month = int(floor(Month));
  *year = int(floor(Year));
  
  *frac = F;

}



void dayofyear(int day, int month, int year, int *doy){
 int lyear = 0;
 int ind;
 
float  daysinmonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};


 if(fmod((float)year,4)==0){
   lyear = 1;
   if(fmod((float)year,100) ==0){ 
     lyear = 0;
     if(fmod((float)year,400)==0){
       lyear=1;
       }
   }
 }
 if(lyear) daysinmonth[1] = 29; /* printf("Year: %d is a leap year!",year); */

 *doy=0;

 for(ind=0;ind<month-1;ind++){
   *doy+= int(daysinmonth[ind]);
 }
 *doy+=day;

}


void geod2geoc(double longitude,  double latitude, double height, const char RefEllips[], double Geoc[], double GeocCart[]){
/* %---------------------------------------------------------------------- */
/* % geod2geoc function        Convert Geodetic coordinates to */
/* %                         Geocentric coordinates, using specified */
/* %                         reference ellipsoid. */
/* % Input  : - Geodetic coordinates. This is three column matrix of the */
/* %            type: [Longitude, Latitude, Height], in case that two */
/* %            column matrix is given, than the height is taken as zero. */
/* %            Units: Longitude and latitude measured in radians while */
/* %            the height is measured in meters above reference ellipsoid. */
/* %          - Reference ellipsoid: */
/* %            'Merit1983' - Merit 1983           a=6378137  1/f=298.257 */
/* %            'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222 */
/* %            'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167 */
/* %            'IAU1976'   - IAU 1976               6378140      298.257 */
/* %            'IAU1964'   - IAU 1964               6378160      298.25 */
/* %            'WGS84'     - WGS 1984 (default) */
/* % Output : - Geocentric coordinates matrix of the type: */
/* %            [longitude, latitude, radius] */
/* %            where longitude and latitude are measured in radians */
/* %            and radius measured in meters from the reference ellipsoid center. */
/* %          - Geocentric cartesian coordinates [x, y, z] in meters. */
/* % Reference : Astronomical Almnach */
/* % Tested : Matlab 5.3 */
/* %     By : Eran O. Ofek */
/* %    URL : http://wise-obs.tau.ac.il/~eran/matlab.html */
/* % SeeAlso: refellipsoid.m, geoc2geod.m */
/* %---------------------------------------------------------------------- */

  double GeodLon,GeodLat,GeodH,C,S,X,Y,Z,R2,GeocLat,GeocRad;
  double A,F;

refellipsoid(RefEllips,&A,&F);

 /* Geodetic position */
 GeodLon = longitude;
GeodLat =  latitude;
GeodH   = height;

C       = pow((pow(cos(GeodLat),2) + pow(((1-F)*sin(GeodLat)),2)),-0.5);
S       = C*pow((1 - F),2);

X       = (A*C + GeodH)*cos(GeodLat)*cos(GeodLon);
Y       = (A*C + GeodH)*cos(GeodLat)*sin(GeodLon);
Z       = (A*S + GeodH)*sin(GeodLat);

R2 = X*X + Y*Y;     /* =(GeocRad*cos(GeoLat)).^2; */
GeocLat = atan(Z/sqrt(R2));
GeocRad = sqrt(R2)/cos(GeocLat);

 Geoc[0] = GeodLon; 
 Geoc[1] = GeocLat;
 Geoc[2] = GeocRad;

 GeocCart[0] = X;
 GeocCart[1] = Y;
 GeocCart[2] = Z;

}


void refellipsoid(const char RefEllips[], double *A, double *F){

/* %---------------------------------------------------------------------- */
/* % refellipsoid function       Return data on a given reference */
/* %                           ellipsoid for the Earth. */
/* % Input  : - Reference ellipsoid: */
/* %            'Merit1983' - Merit 1983           a=6378137  1/f=298.257 */
/* %            'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222 */
/* %            'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167 */
/* %            'IAU1976'   - IAU 1976               6378140      298.257 */
/* %            'IAU1964'   - IAU 1964               6378160      298.25 */
/* %            'WGS84'     - WGS 1984 (default)     6378137      298.257223563 */
/* %            'IERS1989'  - IERS 1989              6378136      298.257 */
/* % Output : - Data matrix: */
/* %            [Equatorial radius (meters), */
/* %             Flattening factor]. */
/* % Reference : Astronomical Almnach */
/* % Tested : Matlab 5.3 */
/* %     By : Eran O. Ofek         June 2000 */
/* %    URL : http://wise-obs.tau.ac.il/~eran/matlab.html */
/* % SeeAlso: geod2geoc.m, geoc2geod.m */
/* %---------------------------------------------------------------------- */

   int ellipsoidchoice(const char RefEllips[]); 

  switch(ellipsoidchoice(RefEllips)){
  case 0:
    /* WGS 1984 */
    *A = 6378137;
    *F = 1/298.257223563;
    break;
  case 1:
    /* IERS 1989 */
    *A = 6378136;
    *F = 1/298.257;
    break;
  case 2:
    /* Merit1983 - Merit 1983           a=6378137  1/f=298.257 */
    *A = 6378137;
    *F = 1/298.257;
    break;
  case 3:
    /* 'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222 */
    *A = 6378137;
    *F = 1/298.257222;
    break;
  case 4:
    /* %'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167 */
    *A = 6378160;
    *F = 1/298.247167;
    break;
  case 5:
    /* %'I*AU1976'   - IAU 1976               6378140      298.257 */
    *A = 6378140;
    *F = 1/298.257;
    break;
  case 6:
    /* %'IAU1964'   - IAU 1964               6378160      298.25 */
    *A = 6378160;
    *F = 1/298.25;
    break;
  default:
    printf("Choice of reference ellipsoid unknown, assuming WGS 1984\n");
    *A = 6378137;
    *F = 1/298.257223563;
  }
 

}
 int ellipsoidchoice(const char RefEllips[]){

 if(strstr(RefEllips,"WGS84")!=NULL) return(0);
 if(strstr(RefEllips,"IERS1989")!=NULL) return(1);
 if(strstr(RefEllips,"Merit1983")!=NULL) return(2);
 if(strstr(RefEllips,"GRS80")!=NULL) return(3);
 if(strstr(RefEllips,"GRS67")!=NULL) return(4);
 if(strstr(RefEllips,"IAU1976")!=NULL) return(5);
 if(strstr(RefEllips,"IAU1964")!=NULL) return(6);

 return(-1);

 }


void mooncool(double JD , double EarthPos[] ,char Algo , double *RA, double *Dec, double *HP){

/* %------------------------------------------------------------------------ */
/* % mooncool function  Calculate low-accuracy Moon Topocentric Equatorial */
/* %                  coordinates (Equinox of date). */
/* % input  : - matrix od dates */
/* %            JD per line. In TT time scale. */
/* %          - [East_Long, North_Lat] of observer in radians. */
/* %            If NaN then calculate geocentric position. */
/* %          - Algorithm: */
/* %            'l' : very low accuracy (default). */
/* %                  0.3 deg in pos. (apparent coordinates). */
/* %                  0.003 deg. in horizontal parallax. */
/* %            'b' : low accuracy ~1' in position. */
/* % output : - vector of RA, in radians. */
/* %          - vector of Dec. in radians. */
/* %          - Vector of horizontal parallax. */
/* %            r = 1/sin(HP)  SD = 0.2725.*HP */
/* %    By  Eran O. Ofek           August 1999 */
/* %------------------------------------------------------------------------ */

  double T,n,SA1,SA2,SA3,SA4,SA5,SA6,Lam;
  double BA1,BA2,BA3,BA4,Bet;
  double CA1,CA2,CA3,CA4,r,l,m,x,y,z;
  double LST,Lat,xt,yt,zt;

T   = (JD - 2451545.0)/36525.0;

 switch(Algo){
 case 'l':

    n  = JD - 2451545.0;
    if (abs(n)>50*365){
    printf("This formulae give good results only between 1950-2050\n");
    }

    SA1 = sin((134.9 + 477198.85*T)/(TO_DEG));
    SA2 = sin((259.2 - 413335.38*T)/(TO_DEG));
    SA3 = sin((235.7 + 890534.23*T)/(TO_DEG));
    SA4 = sin((269.9 + 954397.70*T)/(TO_DEG));
    SA5 = sin((357.5 +  35999.05*T)/(TO_DEG));
    SA6 = sin((186.6 + 966404.05*T)/(TO_DEG));
    Lam = (218.32 + 481267.883*T + 6.29*SA1 - 1.27*SA2 + 0.66*SA3 + 0.21*SA4 - 0.19*SA5 - 0.11*SA6)/(TO_DEG);
    
    BA1 = sin((93.3  + 483202.03*T)/(TO_DEG));
    BA2 = sin((228.2 + 960400.87*T)/(TO_DEG));
    BA3 = sin((318.3 +   6003.18*T)/(TO_DEG));
    BA4 = sin((217.6 - 407332.20*T)/(TO_DEG));
    Bet = (5.13*BA1 + 0.28*BA2 - 0.28*BA3 - 0.17*BA4)/(TO_DEG);
    
    CA1 = cos((134.9 + 477198.85*T)/(TO_DEG));
    CA2 = cos((259.2 - 413335.38*T)/(TO_DEG));
    CA3 = cos((235.7 + 890534.23*T)/(TO_DEG));
    CA4 = cos((269.9 + 954397.70*T)/(TO_DEG));
    *HP  = (0.9508 + 0.0518*CA1 + 0.0095*CA2 + 0.0078*CA3 + 0.0028*CA4)/(TO_DEG);    
    r = 1.0/sin(*HP);
    
    l = cos(Bet)*cos(Lam);
    m = 0.9175*cos(Bet)*sin(Lam) - 0.3978*sin(Bet);
    n = 0.3978*cos(Bet)*sin(Lam) + 0.9175*sin(Bet);
    
    x = r*l;
    y = r*m;
    z = r*n;

    break;


 default:
   printf("Unknown algorithm\n");
   exit(1);

  }

   lst(JD,EarthPos[0],&LST);
   Lat = EarthPos[1];

   xt = x - cos(Lat)*cos(LST*2*PI);
   yt = y - cos(Lat)*sin(LST*2*PI);
   zt = z - sin(Lat);

*RA  = atan2(yt,xt);
*Dec = asin(zt/sqrt(xt*xt + yt*yt + zt*zt));

}




double distsp(double D1, double D2, double R1, double R2){

  double elon;

  elon = acos(sin(D1)*sin(D2) + cos(D1)*cos(D2)*cos(R1-R2));

  return(elon);

}

/*void moon_sky_brightness_old(double JD, double ObjCooRA , double ObjCooDec, double longitude, double latitude, double height, double C_Ext, double Vsky,  double *DeltaV, double *D, double *ObjMoonDist, double *K){

// %----------------------------------------------------------------------- 
// % moon_sky_brightness function     Calculate sky brightness due to the 
// %                                moon for a given date and sky position, 
// %                                by taking into acount the distance from 
// %                                the moon and zenith distance (V mag). 
// % Input  : - Date JD. 
// %          - Object apparent equatorial 
// %            coordinates [RA, Dec] in radians. 
// %          - Observer geodetic position [East_Long, Lat, Height], 
// %            radians and meters above ref ellips. 
// %          - Extinction coef. in V. (default is 0.3mag/airmass). 
// %          - Sky brightness in V. (default is 21.7 mag/sq. arcsec.). 
// % Output : - The change in the V-band sky brightness caused by moonlight. 
// %          - Moon elongation, radians. 
// %          - Object-Moon distance, radians. 
// %          - Moon illuminated fraction. 
// % Reference : Krisciunas, K. and Schaefer, B. 1991 PASP 103, 1033. 

  double Geoc[3];
  double GeocCart[3];
  double MoonRA,MoonDec,MoonHP,MoonR;
  double SunRA,SunDec,SunR,SunSL,SunEquationTime;
  double Alpha,I,Z,Z_Moon;
  double ObjHor_0, ObjHor_1, MoonHor_0, MoonHor_1;
  double I_star, F_Rho, Xz, XzMoon,Bmoon, VskyNL;

// % Geodetic to Geocentric 
// [GeocPos]=geod2geoc(GeodPos,'WGS84'); 

  geod2geoc(longitude,latitude,height, "WGS84", Geoc, GeocCart);

// % Moon & Sun Position 
// [MoonRA,MoonDec,MoonHP] = mooncool(JD,GeocPos); 

mooncool(JD , Geoc ,'l', &MoonRA, &MoonDec, &MoonHP);

// [SunRA,SunDec,SunR]     = suncoo(JD,'a'); 

suncoo(JD, 'a', &SunRA, &SunDec, &SunR, &SunSL, &SunEquationTime);

MoonR = asin(MoonHP)*6378.137/149597870.0; // %AU 

// % Moon Elongation 
*D = distsp(MoonDec,SunDec,MoonRA,SunRA);
Alpha = (TO_DEG)*(PI - *D);


// % Selenographic elongation of the Earth from the Sun 
I = atan2((SunR*sin(*D)),(MoonR - SunR*cos(*D)));

// % Illuminated fraction 
*K = 0.5*(1 + cos(I));


// % convert obj coo. to horiz coo. 
// ObjHor = horiz_coo(ObjCoo,JD,GeodPos,'h'); 

 eq2horiz(ObjCooRA, ObjCooDec, JD, latitude, longitude, 'h', &ObjHor_0, &ObjHor_1);

Z = PI/2 - ObjHor_1;

// % Moon horizontal coo. 
 // MoonHor = horiz_coo([MoonRA,MoonDec],JD,GeodPos,'h');  
 eq2horiz(MoonRA, MoonDec, JD, latitude, longitude, 'h', &MoonHor_0, &MoonHor_1);


Z_Moon = PI/2 - MoonHor_1;


// % Object-Moon distance 
 
// ObjMoonDist = distsp(ObjHor(2),MoonHor(2),ObjHor(1),MoonHor(1)); 
 *ObjMoonDist = distsp(ObjHor_1, MoonHor_1, ObjHor_0, MoonHor_0);


// % moon ilumination 
// I_star = 10.^(-0.4.*(3.84 + 0.026.*abs(Alpha) + 4e-9.*Alpha.^4)); 

    I_star = pow(10,(-0.4*(3.84 + 0.026*fabs(Alpha) + 4e-9*pow(Alpha,4))));

// F_Rho  = (10.^5.36).*(1.06 + cos(ObjMoonDist).^2) + 10.^(6.15 - (TO_DEG).*ObjMoonDist./40); 

 F_Rho  = pow(10,5.36)*(1.06 + pow(cos(*ObjMoonDist),2)) + pow(10,(6.15 - (TO_DEG)* *ObjMoonDist/40));



  Xz     = pow((1 - 0.96*sin(Z)*sin(Z)),(-0.5));

  XzMoon = pow((1 - 0.96*sin(Z_Moon)*sin(Z_Moon)),(-0.5));  //WARNING! ORIGINAL CODE HAS ERROR!  

  // XzMoon = pow((1 - 0.96*sin(Z)*sin(Z)),(-0.5));  THIS LINE SHOULD REPRODUCE ORIGINAL CODE VAULES (WRONG!) 

// % moon sky brigtness in nanoLamberts 
// Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*XzMoon); 

 Bmoon  = F_Rho*I_star*(1 - pow(10,(-0.4*C_Ext*Xz)))*pow(10,(-0.4*C_Ext*XzMoon));


 if(MoonHor_1<0) Bmoon = 0;

// % convert nanLamberts to mag/sq. arcsec. 

// %Bsky   = -(log(Bmoon./34.08) - 20.7233)./0.92104 



// % convert sky brightness to nanoLamberts 

VskyNL = 34.08*exp(20.7233 - 0.92104*Vsky);

*DeltaV = -2.5*log10((Bmoon + VskyNL)/VskyNL);



}*/






























void moon_sky_brightness(double JD, double ObjCooRA , double ObjCooDec, double longitude, double latitude, double height, double C_Ext, double Vsky,  double *DeltaV, double *D, double *ObjMoonDist, double *K){

/* %----------------------------------------------------------------------- */
/* % moon_sky_brightness function     Calculate sky brightness due to the */
/* %                                moon for a given date and sky position, */
/* %                                by taking into acount the distance from */
/* %                                the moon and zenith distance (V mag). */
/* % Input  : - Date JD. */
/* %          - Object apparent equatorial */
/* %            coordinates [RA, Dec] in radians. */
/* %          - Observer geodetic position [East_Long, Lat, Height], */
/* %            radians and meters above ref ellips. */
/* %          - Extinction coef. in V. (default is 0.3mag/airmass). */
/* %          - Sky brightness in V. (default is 21.7 mag/sq. arcsec.). */
/* % Output : - The change in the V-band sky brightness caused by moonlight. */
/* %          - Moon elongation, radians. */
/* %          - Object-Moon distance, radians. */
/* %          - Moon illuminated fraction. */
/* % Reference : Krisciunas, K. and Schaefer, B. 1991 PASP 103, 1033. */

  double Geoc[3];
  double GeocCart[3];
  double MoonRA,MoonDec,MoonHP,MoonR;
  double SunRA,SunDec,SunR,SunSL,SunEquationTime;
  double Alpha,I,Z,Z_Moon;
  double ObjHor_0, ObjHor_1, MoonHor_0, MoonHor_1;
  double I_star, F_Rho, Xz, XzMoon,Bmoon, Bsky, VskyNL;

/* % Geodetic to Geocentric */
/* [GeocPos]=geod2geoc(GeodPos,'WGS84'); */

  geod2geoc(longitude,latitude,height, "WGS84", Geoc, GeocCart);

/* % Moon & Sun Position */
/* [MoonRA,MoonDec,MoonHP] = mooncool(JD,GeocPos); */

mooncool(JD , Geoc ,'l', &MoonRA, &MoonDec, &MoonHP);

/* [SunRA,SunDec,SunR]     = suncoo(JD,'a'); */

suncoo(JD, 'a', &SunRA, &SunDec, &SunR, &SunSL, &SunEquationTime);

MoonR = asin(MoonHP)*6378.137/149597870.0; /* %AU */

/* % Moon Elongation */
*D = distsp(MoonDec,SunDec,MoonRA,SunRA);
Alpha = (TO_DEG)*(PI - *D);


/* % Selenographic elongation of the Earth from the Sun */
I = atan2((SunR*sin(*D)),(MoonR - SunR*cos(*D)));

/* % Illuminated fraction */
*K = 0.5*(1 + cos(I));


/* % convert obj coo. to horiz coo. */
/* ObjHor = horiz_coo(ObjCoo,JD,GeodPos,'h'); */

 eq2horiz(ObjCooRA, ObjCooDec, JD, latitude, longitude, 'h', &ObjHor_0, &ObjHor_1);

Z = PI/2 - ObjHor_1;

/* % Moon horizontal coo. */
 /* MoonHor = horiz_coo([MoonRA,MoonDec],JD,GeodPos,'h');  */
 eq2horiz(MoonRA, MoonDec, JD, latitude, longitude, 'h', &MoonHor_0, &MoonHor_1);


Z_Moon = PI/2 - MoonHor_1;


/* % Object-Moon distance */
 
/* ObjMoonDist = distsp(ObjHor(2),MoonHor(2),ObjHor(1),MoonHor(1)); */
 *ObjMoonDist = distsp(ObjHor_1, MoonHor_1, ObjHor_0, MoonHor_0);


/* % moon ilumination */
/* I_star = 10.^(-0.4.*(3.84 + 0.026.*abs(Alpha) + 4e-9.*Alpha.^4)); */

    I_star = pow(10,(-0.4*(3.84 + 0.026*fabs(Alpha) + 4e-9*hCube(Alpha))));

/* F_Rho  = (10.^5.36).*(1.06 + cos(ObjMoonDist).^2) + 10.^(6.15 - (TO_DEG).*ObjMoonDist./40); */
//10^5.26 = 229086.765276777
 F_Rho  = 229086.765276777*(1.06 + sqr(cos(*ObjMoonDist))) + pow(10,(6.15 - (TO_DEG)* *ObjMoonDist/40.0));



 Xz     = 1.0/sqrt(1 - 0.96*sqr(sin(Z)));

 XzMoon = 1.0/sqrt(1 - 0.96*sqr(sin(Z_Moon)));  /*WARNING! ORIGINAL CODE HAS ERROR!  */

  /* XzMoon = pow((1 - 0.96*sin(Z)*sin(Z)),(-0.5));  THIS LINE SHOULD REPRODUCE ORIGINAL CODE VAULES (WRONG!) */

/* % moon sky brigtness in nanoLamberts */
/* Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*XzMoon); */

 Bmoon  = F_Rho*I_star * (1 - pow(10,(-0.4*C_Ext*Xz))) * pow(10,(-0.4*C_Ext*XzMoon));


 if(MoonHor_1<0) Bmoon = 0;

/* % convert nanLamberts to mag/sq. arcsec. */

/* %Bsky   = -(log(Bmoon./34.08) - 20.7233)./0.92104 */


//Night sky brightness (no moon - airglow only)

//Assume this is the background specified in the detector file (at zenith)

/* % convert sky brightness to nanoLamberts */
VskyNL = 34.08*exp(20.7233 - 0.92104*Vsky);
Bsky = VskyNL * pow(10,-0.4*C_Ext*(Xz-1)) * Xz;
 




*DeltaV = -2.5*log10((Bmoon + Bsky)/VskyNL);



}






















void mag2flux(double mag, char band, double lambda0_def, double deltalambda_def, double f_nu0_def,double *flux){
  double lambda0,deltalambda,f_nu0;
  double frac, logfnu0;

/* Convert magnitude <mag> in band <band> to photon flux value <flux>. */
/* If <band>='X', use defined lambda0_def, deltalambda_def and f_nu0_def */
/* values. */
  
/* Standard passband parameters from Zombeck */

/* Added a Y-band with lambda, delta lambda taken from EUCLID SYS_PERF_REF 
   iss6 document, and value of f_nu0 guessed */

  /* Pre-computed: frac = deltalambda/lambda0; logfnu0 = log10(f_nu0) */

/* Compare to http://www.astro.soton.ac.uk/~rih/applets/MagCalc.html */

switch(band){
 case 'U':
   lambda0 = 0.365; deltalambda= 0.068; f_nu0 = 1.90e-23;
   frac = 0.186301369863014; logfnu0 = -22.721246399047168;
   break;
 case 'V':
lambda0 = 0.44; deltalambda= 0.098; f_nu0 = 4.27e-23;
 frac = 0.222727272727273; logfnu0 =-22.369572124974972;
   break;
 case 'B':
lambda0 = 0.55; deltalambda= 0.089; f_nu0 = 3.67e-23;
 frac =0.161818181818182; logfnu0 =-22.435333935747913;
   break;
 case 'R':
lambda0 = 0.70; deltalambda= 0.22; f_nu0 = 2.84e-23;
 frac = 0.314285714285714; logfnu0 = -22.546681659952966;
   break;
 case 'I':
lambda0 = 0.90; deltalambda= 0.24; f_nu0 = 2.25e-23;
 frac = 0.266666666666667; logfnu0 = -22.647817481888637;
   break;
 case 'Y':
   lambda0 = 1.033; deltalambda=0.226; f_nu0 = 2.0e-23;
 frac = 0.218780252; logfnu0 = -22.698970004;
   break;
 case 'J':
lambda0 = 1.25; deltalambda= 0.3; f_nu0 = 1.65e-23;
 frac = 0.240; logfnu0 = -22.782516055786090;
   break;
 case 'H':
lambda0 = 1.65; deltalambda= 0.4; f_nu0 = 1.07e-23;
 frac = 0.242424242424242; logfnu0 = -22.970616222314789;
   break;
 case 'K':
lambda0 = 2.2; deltalambda= 0.6; f_nu0 = 6.73e-24;
 frac = 0.272727272727273; logfnu0 = -23.171984935776024;
   break;
 case 'L':
lambda0 = 3.6; deltalambda= 1.2; f_nu0 = 2.69e-24;
 frac = 0.333333333333333; logfnu0 = -23.570247719997596;
   break;
 case 'X':
lambda0 = lambda0_def; deltalambda= 0.068; f_nu0 = f_nu0_def;
 frac = deltalambda/lambda0_def; logfnu0 = log10(f_nu0_def);
   break;

 default:
   printf("(astroFns:mag2flux) Unknown passband: %c .\n",band); exit(1);
 }

 *flux = pow(10,-0.4*mag + logfnu0)* 1.51e33 * frac; /* Photons per second per metre^2 */


 /* HARD CODED - APPLIES TO EUCLID ONLY!!!!!!! */
 
 switch(band)
   {
   case 'H':
     frac = 121.122;
     break;
   case 'Y':
     frac = 79.1129;
     break;
   case 'J':
     frac = 65.4916;
     break;
   case 'R':
   case 'I':
     frac = 49.3396;
     break;
   default:
     return;
     
   }
 *flux = frac*pow(10,-0.4*(mag-20.0));

}

int bandcode(char band)
{
  int colidx;

  switch(band)
    {
    case 'R':
      colidx=0;
      break;
    case 'I':
      colidx=1;
      break;
    case 'Y':
      colidx=2;
      break;
    case 'J':
      colidx=3;
      break;
    case 'H':
      colidx=4;
      break;
    default:
      colidx=4;
      break;
    }

  return colidx;
}

//in degrees
void gal2eclip(double l, double b, double* lambda, double* beta)
{
  const double sinbetaNGP = 0.497125407;
  const double cosbetaNGP = 0.867678702;
  const double l1=6.38;
  const double lambda0=270.02;
 
  double sinb=sin(b*TO_RAD);
  double cosb=cos(b*TO_RAD);
  double sinlml1=sin((l-l1)*TO_RAD);
  double coslml1=cos((l-l1)*TO_RAD);
 
  double coslml0, sinlml0, sinbeta, cosbeta;
 
 
  sinbeta = sinb*sinbetaNGP + cosb*cosbetaNGP*sinlml1;
  cosbeta = sqrt(1.0-sinbeta*sinbeta);
  coslml0 = coslml1*cosb/cosbeta;
  sinlml0 = (-sinb*cosbetaNGP + cosb*sinbetaNGP*sinlml1)/cosbeta;
 
  *lambda = atan2(sinlml0,coslml0)*TO_DEG + lambda0;
  while(*lambda<0){
    *lambda += 360;
  }
  while(*lambda>360){
    *lambda -= 360;
  }
  *beta = asin(sinbeta)*TO_DEG;
   
}
 
//in degrees
void eq2eclip(double ra, double dec, double* lambda, double* beta)
{
  const double sineps = 0.397772494;
  const double coseps = 0.917484083;
 
  double sindec = sin(dec*TO_RAD);
  double cosdec = cos(dec*TO_RAD);
  double sinra = sin(ra*TO_RAD);
  double cosra = cos(ra*TO_RAD);
 
  double cosl, sinl, sinbeta, cosbeta;
 
  sinbeta = sindec*coseps - cosdec*sineps*sinra;
  cosbeta = sqrt(1.0-sinbeta*sinbeta);
  cosl = cosra*cosdec/cosbeta;
  sinl = (sindec*sineps + cosdec*coseps*sinra)/cosbeta;
 
  *lambda = atan2(sinl,cosl)*TO_DEG;
  while(*lambda<0){
    *lambda += 360;
  }
  while(*lambda>360){
    *lambda -= 360;
  }
  *beta = asin(sinbeta)*TO_DEG;
   
}
 
//in degrees
double solarlambda(double jdate)
{
  double ra, dec, R, SL, EqTime;
  double lambda, beta;
 
  suncoo(jdate, 'j', &ra, &dec, &R, &SL, &EqTime);
   
  eq2eclip(ra*TO_DEG,dec*TO_DEG,&lambda,&beta);
 
  return lambda;
   
}
