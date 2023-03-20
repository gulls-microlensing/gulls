/*! \file 
\brief Functions related to determining blending

These functions are used to estimate the blending fraction f_S
for the event. 
 */

#include "findBlending.h"
#include "blendHist_s1d1.h"
#include "blendHist_s1d3.h"
#include "blendHist_s1d5.h"
#include "blendHist_s4d1.h"
#include "blendHist_s4d3.h"
#include "blendHist_s4d5.h"
#include "blendHist_s7d1.h"
#include "blendHist_s7d3.h"
#include "blendHist_s7d5.h"

#define DEBUGVAR 0

/*! Read in the blending histograms from header files: 
  blendHist_sXdY.h 

 where X in [1,4,7] is seeing and Y in [1 3 5] is  
 stellar density.  */
void initBlending(struct hist (*BlendData)[3]){
  char str[200];
 
  blendHist_s1d1(BlendData);
  blendHist_s1d3(BlendData);
  blendHist_s1d5(BlendData);
  blendHist_s4d1(BlendData);
  blendHist_s4d3(BlendData);
  blendHist_s4d5(BlendData);
  blendHist_s7d1(BlendData);
  blendHist_s7d3(BlendData);
  blendHist_s7d5(BlendData);


  sprintf(str,"%s","Blending histograms"); 
  fmtline(str,WIDTH,"PARSED"); 
 
}


void findBlending(struct obsfilekeywords World[], int numobservatories, struct event *Event,long *idum,struct hist (*BlendData)[3], struct filekeywords* Paramfile)
{
  void stellardensity(double longitude, double latitude, int *stellardensityidx);
  void computefs (double dice, int obsidx, int colidx, int stellardensityidx, struct obsfilekeywords World[], struct event *Event,long *idum,struct hist (*BlendData)[3]);

  int stellardensityidx,obsidx,colidx,filter;
  double dice;
  double fs;

  double Fsource, Flens, Fblend;

  dice = ran2(idum);
  
  stellardensity(Event->l, Event->b, &stellardensityidx);  /*Find stellar density in direction of event */

  switch(Paramfile->primaryColour[0])
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

  /* Use the same fs for all bands and observatories */
  computefs(dice, 0, colidx, stellardensityidx, World, Event,idum,BlendData);

  fs = Event->fs[0];

  for(obsidx=0;obsidx<numobservatories;obsidx++)
    { 
      /* This version has correlated but different blending in each band
	 computefs(dice, obsidx, colidx, stellardensityidx, World, Event,idum,BlendData); 
	 fs = Event->fs[0];
      */

      /* Include the effect of the lens */
      filter = World[obsidx].filter;

      Fsource = pow(10,-0.4*(Event->m_s.m[filter]-20));
      Flens   = pow(10,-0.4*(Event->m_l.m[filter]-20));
      Fblend  = Fsource*(1-fs)/fs;

      Event->fs[obsidx] = Fsource/(Fsource+Flens+Fblend);
	  
    } 

}

/*! This is where we figure out the stellar density for the  
 event at Galactic l and Galactic b and return a value in  
 stellardensity which is an index to a table of cumulative 
 frequencies. 

 CURRENTLY HARD-CODED TO RETURN STELLARDENSITYIDX CORRESPONDING TO MEDIUM
 DENSITY FIELDS

*/
void stellardensity(double longitude, double latitude, int *stellardensityidx)
{
  /*Hard coded to return low-stellar density field results*/
  *stellardensityidx = 1;
}


/*! This is where we would look up the cumulative frequency table of fs,  
  given the stellar density (indexed by stellardensityidx) and each 
  site's mean seeing */
void computefs (double dice, int obsidx, int colidx, int stellardensityidx, struct obsfilekeywords World[], struct event *Event,long *idum, struct hist (*BlendData)[3])
{
  double a,b,c;
  int seeingidx=9;
  int mdx=100;
  int knd;
  int idx=0;
  double interpval,x1,x2,y1,y2;

  if(DEBUGVAR) printf("blending profile = %d\n", World[obsidx].blendingProfile);

  switch(World[obsidx].blendingProfile)
    {
      
    case 0:    /* No blending */
      if(DEBUGVAR) printf("Space observatory. Setting blending to 1.0\n");
      /*fprintf(logfile_ptr,"BLENDING: Space observatory - setting blending to 1.0\n");*/
      Event->fs[obsidx]=1.0;
      return;

    case 11:
      seeingidx = 0;
      stellardensityidx = 0;
      break;

    case 13:
      seeingidx = 0;
      stellardensityidx = 1;
      break;

    case 15:
      seeingidx = 0;
      stellardensityidx = 2;
      break;

    case 41:
      seeingidx = 1;
      stellardensityidx = 0;
      break;

    case 43:
      seeingidx = 1;
      stellardensityidx = 1;
      break;

    case 45:
      seeingidx = 1;
      stellardensityidx = 2;
      break;

    case 71:
      seeingidx = 2;
      stellardensityidx = 0;
      break;

    case 73:
      seeingidx = 2;
      stellardensityidx = 1;
      break;

    case 75:
      seeingidx = 2;
      stellardensityidx = 2;
      break;

    case 1:
      stellardensityidx = 0;
    case 3:
      stellardensityidx = 1;
    case 5:
      stellardensityidx = 2;
    default:

      if(World[obsidx].space)
	{
	  if(DEBUGVAR) printf("Space observatory. Setting blending to 1.0\n");
	  Event->fs[obsidx] = 1.0;
	  /*fprintf(logfile_ptr,"BLENDING: Space observatory - setting blending to 1.0\n");*/
	  return;
	}

      a = fabs(World[obsidx].meanseeing-0.7);
      b = fabs(World[obsidx].meanseeing-1.05);
      c = fabs(World[obsidx].meanseeing-2.10);

      if(a<b && a<c) seeingidx = 0;
      if(b<a && b<c) seeingidx = 1;
      if(c<b && c<a) seeingidx = 2;

    };

  if(DEBUGVAR) printf("seeingidx: %d  stellardensityidx: %d\n",seeingidx, stellardensityidx);

  /*find the event magnitude range*/
  for(knd=0;knd<NUM_RANGES_BLENDHIST;knd++)
    {
      switch(knd)
	{
	case 0:
	  if(Event->m_s.m[colidx] <
	     BlendData[seeingidx][stellardensityidx].rangeMax[knd])
	    mdx=knd;
	  break;

	case NUM_RANGES_BLENDHIST-1:
	  if(mdx==100) mdx=knd;
	  break;
	  
	default:
	  if(Event->m_s.m[colidx] >= 
	     BlendData[seeingidx][stellardensityidx].rangeMin[knd] 
	     && Event->m_s.m[colidx] <
	     BlendData[seeingidx][stellardensityidx].rangeMax[knd])
	    mdx=knd;
	}
    }

  idx =0;
  while(dice>BlendData[seeingidx][stellardensityidx].yHist[mdx][idx])
    {
      /*if(DEBUGVAR) 
	printf("%d %f %f\n", idx, dice, 
	BlendData[seeingidx][stellardensityidx].yHist[mdx][idx]);*/

      idx++;
    }
  
  if(idx==0)
    {
      
      Event->fs[obsidx]
	= BlendData[seeingidx][stellardensityidx].xHist[mdx][idx];
      if(DEBUGVAR)
	{
	  printf("Dice roll gives idx=0 Using value at idx=0, i.e. no interpolation\n");
	  printf("BlendData[%d][%d].xHist[%d][%d] = %f\n",seeingidx,stellardensityidx,mdx,idx,BlendData[seeingidx][stellardensityidx].xHist[mdx][idx]);
	}
      return;
    }

  y1 = BlendData[seeingidx][stellardensityidx].yHist[mdx][idx-1];
  y2 = BlendData[seeingidx][stellardensityidx].yHist[mdx][idx];
  x1 = BlendData[seeingidx][stellardensityidx].xHist[mdx][idx-1];
  x2 = BlendData[seeingidx][stellardensityidx].xHist[mdx][idx];
  
  if(y2>99.9) y2 = 1.0;
  
  interpval = (x2-x1)/(y2-y1)*(dice-y1)+x1;
  
  /* The following line would be a function of stellardensityidx and */
  /* seeingidx, and is the direct extraction from the histogram data */
  Event->fs[obsidx] =
    BlendData[seeingidx][stellardensityidx].xHist[mdx][idx];
  
  if(DEBUGVAR) printf("fs direct = %f, mdx = %d, idx = %d\n",Event->fs[obsidx],mdx,idx);

  /* Using interpolated value */
  Event->fs[obsidx] = interpval;
  if(DEBUGVAR) printf("fs interp = %f, mdx = %d, idx = %d\n",Event->fs[obsidx],mdx,idx);
  
  if(DEBUGVAR)
    {
      printf("Observatory: %d Filter: %d Mean seeing: %f. Seeingidx = %d\n", obsidx, colidx, World[obsidx].meanseeing,seeingidx);
      printf("Stellar density index: %d\n",stellardensityidx);
      printf("Source apparant mag: %f. Corresponds to range: %d (%f -- %f)\n", Event->m_s.m[colidx], mdx, BlendData[seeingidx][stellardensityidx].rangeMin[mdx], BlendData[seeingidx][stellardensityidx].rangeMax[mdx]);
      printf("%f <= Dice = %f < %f.   Corresponds to index: %d\n", BlendData[seeingidx][stellardensityidx].yHist[mdx][idx-1], dice, BlendData[seeingidx][stellardensityidx].yHist[mdx][idx], idx);
      
      printf("BlendData[%d][%d].xHist[%d][%d] = %f  BlendData[%d][%d].xHist[%d][%d] = %f\n", seeingidx, stellardensityidx, mdx, idx-1, BlendData[seeingidx][stellardensityidx].xHist[mdx][idx-1], seeingidx, stellardensityidx, mdx, idx, BlendData[seeingidx][stellardensityidx].xHist[mdx][idx]);
      printf("Interpolated value for blending at dice = %f is %f\n", dice, interpval);
      
    }
  
  return;

}
