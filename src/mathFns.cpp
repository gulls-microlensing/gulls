#include "mathFns.h"


float fsgn(double a)         /* returns the sign of a float:             */
                            /* 1 if positive, -1 if negative, 0 if zero */
{
    if (a>0.0)
        return (1.0);
    else
    {
        if (a<0.0)
            return (-1.0);
        else
            return (0.0);
    }
}

int dsgn(double a)         /* returns the sign of a float:             */
                            /* 1 if positive, -1 if negative, 0 if zero */
{
    if (a>0.0)
        return (1);
    else
    {
        if (a<0.0)
            return (-1);
        else
            return (0);
    }
}


double dmin(double a, double b){

  if(a<b){
    return(a);
  }

  if(a==b){
    return(a);
  }
  if(a>b){
    return(b);
  }

  return(0);
}

double dmax(double a, double b){

  if(a<b){
    return(b);
  }

  if(a==b){
    return(b);
  }
  if(a>b){
    return(a);
  }

  return(0);
}

