#ifndef SOURCE_CLD
#define SOURCE_CLD

#include "cd.h"
#include "src_base.h"
#include "constants.h"

///////////////////////////////////////////////////////////////////////
//
//     src_cld
//
//        Source class for a circular limb darkened source
//        
//
///////////////////////////////////////////////////////////////////////

class src_cld: public src_base
{

 private:

 public:
  
  src_cld(double rs_, double ldGamma_)
    {
      rs=rs_;
      ldGamma=ldGamma_;
      
      twoCoeff=0.5*(1.0-0.2*ldGamma);
      fourCoeff=(1.0/3.0)*(1.0-11.0*ldGamma/35.0);
      area=pi*sqr(rs);
    }

  src_cld()
    {
      src_cld(1.0e-3,0.0);
    }

  void setgamma(double ldGamma_)
    {
      ldGamma=ldGamma_;

      twoCoeff=0.5*(1.0-0.2*ldGamma);
      fourCoeff=(1.0/3.0)*(1.0-11.0*ldGamma/35.0);
    }

  ~src_cld(){};

  inline double limbDarkening(cd& zssc)
  {
    return (ldGamma==0.0?1.0:(1.0 - ldGamma*(1.0 - 1.5*sqrt(1.0-sqr(abs(zssc)/rs)))));
  }

  inline bool inSource(cd& zsc,cd& zs)
  {
    return (abs(zs-zsc)<rs ? true : false);
  }

  inline bool inSource(cd& zssc)
  {
    return (abs(zssc)<rs ? true : false);
  }
};


#endif /*SOURCE_CLD*/
