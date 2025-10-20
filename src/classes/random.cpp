#include "random.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <cmath>

// Fallback implementations for numerical recipes random functions
// Uses GSL when Numerical Recipes are not available
// The actual implementations should be provided by the user

// Global GSL random number generator for fallback
static gsl_rng* gsl_rng_fallback = nullptr;

// Initialize fallback RNG if not already done
static void init_fallback_rng() {
    if (gsl_rng_fallback == nullptr) {
        gsl_rng_fallback = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gsl_rng_fallback, 12345); // Default seed
    }
}

// Cleanup function (called at program exit)
static void cleanup_fallback_rng() {
    if (gsl_rng_fallback != nullptr) {
        gsl_rng_free(gsl_rng_fallback);
        gsl_rng_fallback = nullptr;
    }
}

// Register cleanup function
static int dummy = (atexit(cleanup_fallback_rng), 0);

double ran1(long *idum) {
    init_fallback_rng();
    // Only reseed if idum is negative (NR convention for initialization)
    if (idum && *idum < 0) {
        gsl_rng_set(gsl_rng_fallback, -(*idum));
        *idum = 1; // Mark as initialized
    }
    return gsl_rng_uniform(gsl_rng_fallback);
}

double ran2(long *idum) {
    init_fallback_rng();
    // Only reseed if idum is negative (NR convention for initialization)
    if (idum && *idum < 0) {
        gsl_rng_set(gsl_rng_fallback, -(*idum));
        *idum = 1; // Mark as initialized
    }
    return gsl_rng_uniform(gsl_rng_fallback);
}

double ran0(long *idum) {
    return ran2(idum);
}

double gasdev(long *idum) {
    init_fallback_rng();
    // Only reseed if idum is negative (NR convention for initialization)
    if (idum && *idum < 0) {
        gsl_rng_set(gsl_rng_fallback, -(*idum));
        *idum = 1; // Mark as initialized
    }
    return gsl_ran_gaussian(gsl_rng_fallback, 1.0);
}

double gammln(double xx) {
    // Use GSL's log gamma function
    return gsl_sf_lngamma(xx);
}

double gammp(double a, double x) {
    // Use GSL's incomplete gamma function P(a,x)
    return gsl_sf_gamma_inc_P(a, x);
}

double gammq(double a, double x) {
    // Use GSL's incomplete gamma function Q(a,x) = 1 - P(a,x)
    return gsl_sf_gamma_inc_Q(a, x);
}

int randint(int min, int max, long *seed) {
    init_fallback_rng();
    // Only reseed if seed is negative (NR convention for initialization)
    if (seed && *seed < 0) {
        gsl_rng_set(gsl_rng_fallback, -(*seed));
        *seed = 1; // Mark as initialized
    }
    return min + gsl_rng_uniform_int(gsl_rng_fallback, max - min + 1);
}

double poisson(double mean, long *seed) {
    init_fallback_rng();
    // Only reseed if seed is negative (NR convention for initialization)
    if (seed && *seed < 0) {
        gsl_rng_set(gsl_rng_fallback, -(*seed));
        *seed = 1; // Mark as initialized
    }
    return gsl_ran_poisson(gsl_rng_fallback, mean);
}