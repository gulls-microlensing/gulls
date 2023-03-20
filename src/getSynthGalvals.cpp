/*! \file 
\brief Functions to read in Galaxy data

This file contains functions that reads in the synthetic Galaxy model data.
 */

#include "getSynthGalvals.h"
#include "strfns.h"
#include "random.h"
#include "constdefs.h"

#define DEBUGVAR 0

double dl,ds,l,b,vperp,Rs,ml,M_s_I,M_l_I;

/*! Main code to read and parse Galaxy model data*/
void getSynthGalvals(struct event *Event, struct galaxy *Galaxy, int sdx, long *idum)
{

  if(0)
    {
      /* THESE VALUES WILL ULTIMATELY COME FROM THE SYNTHETIC MODEL */
      /* A LOT OF WORK STILL TO BE DONE HERE  */
      
      dl = 6000;               /* Lens and source distance in pc */
      ds = 8000;
      l = 0.67;            /* Galactic longitude and latitude in degrees */
      b = -3.2362;
      Rs = 1.0;            /* Source star size in units of R_sol */
      ml = 0.3;             /* Mass of lens in solar units */
      vperp = 200;           /* lens transverse velocity in km/s */
 
      /* NEEDS TO BE EXPANDED TO INCLUDE MORE PASSBANDS */
      M_s_I = 0.0;         /* Absolute magnitude I band for source */
      M_l_I = 12.0;        /* Absolute magnitude I band for lens */


      /* ATTEMPT TO TEST FOR VARIOUS PARAMETER VALUES */
      Rs = 1.0 * ran2(idum) + 1.0;
      M_s_I += 5*ran2(idum);
      vperp += 50*gasdev(idum);
      dl+= 500*gasdev(idum);
      ds+= 500*gasdev(idum);
  
    }

  if(DEBUGVAR) printf("getSynthGalvals: *Event = %p\n",Event);
  if(DEBUGVAR) printf("getSynthGalvals: *Galaxy = %p\n",Galaxy);
  if(DEBUGVAR) printf("getSynthGalvals: sdx = %d\n",sdx);
  if(DEBUGVAR) printf("getSynthGalvals: idum = %d\n",(int)*idum);

  /* LOAD MODEL VALUES INTO EVENT STRUCTURE */
  /* Event->dl = Galaxy->dL[sdx] * 1000; */
  /* Event->ds = Galaxy->dS[sdx] * 1000; */
  /* Event->l = Galaxy->l[sdx]; */

  /* if(Event->l > 180.0) Event->l -= 360.0; */

  /* Event->b = Galaxy->b[sdx]; */
  /* Event->Rs = Galaxy->Rs[sdx]; */
  /* Event->ml = Galaxy->mL[sdx]; */
  /* Event->vperp = Galaxy->vperp[sdx];/\* * 47.4057581; now not needed *\/ */
  /* Event->w = Galaxy->w[sdx]; */
  
  /* Event->tE = Galaxy->tE[sdx]; */
  
  /* Event->eventArea = Galaxy->eA[sdx]; */

  /* Event->m_l = Galaxy->m_L[sdx]; */
  /* Event->m_s = Galaxy->m_S[sdx]; */

  /* Event->mp = Galaxy->mp[sdx]*MEARTH;  /\* planet mass in solar masses *\/ */
  /* Event->semiMajor = Galaxy->a[sdx]; */
  /* Event->inc = Galaxy->inc[sdx]*PI/180.0; */
  /* Event->phase = Galaxy->phase[sdx]*PI/180.0; */

  /* Event->id = Galaxy->id[sdx]; */
  /* Event->fieldno = Galaxy->fieldno[sdx]; */

}

