#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "singleLens.h"
#include<time.h>
#include<vector>

#include<fstream>
#include <sstream>
#include<iomanip>
#include <sys/stat.h>
#define DEBUGVAR 1
// This assumes Event->vbm has already been initialized and configured

void lightcurveGenerator(struct filekeywords* Paramfile, struct event *Event, struct obsfilekeywords World[], struct slcat *Sources, struct slcat *Lenses, ofstream& logfile_ptr)
{
    char str[100];
    double rs = Event->rs;
    double u0 = Event->u0;
    double alpha = Event->alpha * TO_RAD;
    double q = Event->params[QQ];
    double a = Event->params[SS];
    double m1 = 1.0 / (1.0 + q);
    double cosa = cos(alpha);
    double sina = sin(alpha);
    double VBM_origin = (1.0 - m1) * (-a);
    vector<int> obsoffset(Paramfile->numobservatories,0);
    Event->Amax=-1;
    Event->umin=1e50;
    Event->lcerror=0;
    int errflag=0;
    int obsidx;
    double lim_gamma=Event->gamma;
    if(Paramfile->verbosity>=3)
	    std::cout 
		    << "At lightcurveGenerator start, Tol=" 
		    << Event->vbm->Tol 
		    << ", RelTol=" 
		    << Event->vbm->RelTol 
		    << std::endl;
     
    //if the event is saturated in each band, no need to calculate the lightcurve
    if(Event->nepochs==0 || Event->allsat)
     {
      return;
     }
    // Resize arrays to match number of epochs
    if(Paramfile->verbosity>=4)
    	{
            Event->vbm_rootaccuracy.resize(Event->nepochs);
            Event->vbm_squarecheck.resize(Event->nepochs);
            Event->vbm_therr.resize(Event->nepochs);
    	}
    vector<int> idxshift;
    int shiftedidx;
    for(obsidx=0;obsidx<Paramfile->numobservatories;obsidx++)
       idxshift.push_back(Event->nepochsvec[obsidx]);

    //Calculate the lightcurve
    for (int idx = 0; idx < Event->nepochs; ++idx) {
        obsidx = Event->obsidx[idx];
        shiftedidx=idx-idxshift[obsidx];
        double amp = 0.0;
        if (Paramfile->identicalSequence && obsidx > 0) {
            // Lightcurve is identical from observatory to observatory
            amp = Event->Atrue[shiftedidx];
        } else {

            double tt = (Event->epoch[idx] - Event->t0) / Event->tE_r;
            double uu = u0;

            if (Paramfile->pllxMultiplyer) {
                tt += Event->pllx[obsidx].tshift[shiftedidx];
                uu += Event->pllx[obsidx].ushift[shiftedidx];
            }
            Event->umin=min(Event->umin,qAdd(tt,uu));
            
	    double xsCoM = tt * cosa - uu * sina + VBM_origin;
            double ysCenter = tt * sina + uu * cosa;

            Event->xs[idx] = xsCoM;
            Event->ys[idx] = ysCenter;
            Event->xl1[idx] = VBM_origin;
            Event->yl1[idx] = 0.0;
            Event->xl2[idx] = VBM_origin + a;
            Event->yl2[idx] = 0.0;
	    Event->vbm->a1 = lim_gamma;
	    amp = Event->vbm->BinaryMag2(a, q, xsCoM, ysCenter, rs);
	 }
        Event->Atrue[idx] = amp;

        if (Paramfile->verbosity >= 4) {
             // Save VBM diagnostics
             Event->vbm_rootaccuracy[idx] = Event->vbm->rootaccuracy;
             Event->vbm_squarecheck[idx] = Event->vbm->squarecheck;
             Event->vbm_therr[idx] = Event->vbm->therr;
            }
        

        // Keep track of highest magnification
        if (amp > Event->Amax) {
            Event->Amax = amp;
            Event->peakpoint = idx;
        }
    }
  //if there has been an error - try the backup generator
  if(Event->lcerror)
    {
      Event->lcerror=0;
      backupGenerator(Paramfile, Event, World, Sources, Lenses, logfile_ptr);
    }
}

