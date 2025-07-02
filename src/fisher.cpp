#include <sstream>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_deriv.h>

#include "lightcurveGenerator.h"
#include "backupGenerator.h"
#include "astroFns.h"
#include "VBMicrolensingLibrary.h"
#include "fs.h"
#include "src_cld.h"
#include "lens_binary.h"
#include "cd.h"

#define DEBUGVAR 0

// Context struct for GSL derivative wrapper
struct DerivContext {
    int param;
    int idx;
    filekeywords* Paramfile;
    event* Event;
    obsfilekeywords* World;
};

// GSL derivative wrapper using context struct
static double deriv_wrapper(double x, void* p) {
    DerivContext* ctx = static_cast<DerivContext*>(p);

    filekeywords* Paramfile = ctx->Paramfile;
    event* Event = ctx->Event;
    obsfilekeywords* World = ctx->World;
    int param = ctx->param;
    int idx = ctx->idx;

    double logq = log10(Event->params[QQ]);
    double logs = log10(Event->params[SS]);
    double logtE = log10(Event->tE_r);
    double piEN = Event->piEN;
    double piEE = Event->piEE;

    double t0 = Event->t0, tE = pow(10, logtE), u0 = Event->u0, alpha = Event->alpha * TO_RAD;
    double rs = Event->rs;

    switch (param) {
        case 0: t0 += x; break;
        case 1: logtE += x; tE = pow(10, logtE); break;
        case 2: u0 += x; break;
        case 3: alpha += x; break;
        case 4: logs += x; break;
        case 5: logq += x; break;
        case 6: rs = pow(10, log10(Event->rs) + x); break;
        case 7: piEN += x; break;
        case 8: piEE += x; break;
    }

    double q = pow(10, logq);
    double a = pow(10, logs);
    double m1 = 1.0 / (1.0 + pow(10, logq));
    double cosa = cos(alpha), sina = sin(alpha);
    double origin = (1.0 - m1) * (-a);

    int obsidx = Event->obsidx[idx];
    int shiftedidx = idx - Event->nepochsvec[obsidx];

    if (Paramfile->pllxMultiplyer && (param == 7 || param == 8 || param == 1)) {
        for (int i = 0; i < Paramfile->numobservatories; i++) {
            Event->pllx[i].provide_observables_NE(piEN, piEE, tE);
            Event->pllx[i].compute_tushifts();
        }
    }

    double tt = (Event->epoch[idx] - t0) / tE;
    double uu = u0;
    if (Paramfile->pllxMultiplyer) {
        tt += Event->pllx[obsidx].tshift[shiftedidx];
        uu += Event->pllx[obsidx].ushift[shiftedidx];
    }

    double xs = tt * cosa - uu * sina + origin;
    double ys = tt * sina + uu * cosa;
    double lim_gamma=Paramfile->LD_GAMMA;
    Event->vbm->a1=lim_gamma;

    return Event->vbm->BinaryMag2(a, q, xs, ys, rs);
}

void fisherMatrix(filekeywords* Paramfile, event* Event,
                  obsfilekeywords World[], slcat* Sources, slcat* Lenses) {

    if (DEBUGVAR) {
        std::cout << "Calculating Fisher matrix for event " << Event->id << std::endl;
    }

    const int Nparams = (Paramfile->pllxMultiplyer == 0) ? 7 : 9;
    const int Nobsparams = 2 * Paramfile->numobservatories;
    const int Ntotparams = Nparams + Nobsparams;

    std::vector<double> step(Ntotparams);
    for (int i = 0; i < Nparams; ++i) {
        double base = (i < Event->params.size()) ? std::max(fabs(Event->params[i]), 1.0) : 1.0;
        if (i == 6) step[i] = 1e-2 * base;
        else if (i == 7 || i == 8) step[i] = 1e-1 * base + 1e-5;
        else if (i == 3) { step[i] = (1e-4 * base) * TO_RAD; } // alpha is stored in degrees; multiply by TO_RAD so δalpha is in radians 
	else step[i] = 1e-4 * base;
    }

    Event->dF.clear();
    Event->dF.resize(Ntotparams * Event->nepochs, 0.0);

    for (int param = 0; param < Nparams; ++param) {
        DerivContext ctx{param, 0, Paramfile, Event, World};
        gsl_function F;
        F.function = &deriv_wrapper;
        F.params   = &ctx;

	for (int idx = 0; idx < Event->nepochs; ++idx) {
		ctx.idx = idx;
		double result = 0.0, abserr = 0.0;
		gsl_deriv_central(&F, 0.0, step[param], &result, &abserr);
		int obsidx = Event->obsidx[idx];
		double fs = Event->fs[obsidx];
		double snr = fabs(result) / (abserr + 1e-12);
		bool suppress = std::isnan(result) || std::isnan(abserr) ||
                            ((param == 6) ? (snr < 1.0) : (snr < 3.0));
		bool in_anomaly_window = fabs(Event->epoch[idx] - Event->t0) < 0.2 * Event->tE_r; //should we keep this for Earth mass planets?
		int index = param * Event->nepochs + idx;
		if (suppress && !in_anomaly_window) {
			double interp = 0.0;
			int count = 0;
			for (int offset = 1; offset <= 3; ++offset) {
				for (int j = -1; j <= 1; j += 2) {
					int neighbor = idx + j * offset;
					if (neighbor >= 0 && neighbor < Event->nepochs) {
						double val = Event->dF[param * Event->nepochs + neighbor];
						if (!std::isnan(val)) {
							interp += val;
							++count;
						}
					}
				}
			}
			Event->dF[index] = (count > 0) ? (interp / count) : result * fs;
		} else {
			Event->dF[index] = result * fs;
		}
	}
    }

    for (int idx = 0; idx < Event->nepochs; ++idx) {
        int obsidx = Event->obsidx[idx];
        int base = (Nparams + obsidx * 2) * Event->nepochs + idx;
        Event->dF[base] = Event->Atrue[idx];
        Event->dF[base + Event->nepochs] = (Event->Atrue[idx] - 1.0) / Event->fs[obsidx];
    }

}

