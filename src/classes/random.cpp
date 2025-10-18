#include "random.h"
#include <cstdlib>
#include <cmath>

// Stub implementations for numerical recipes random functions
// These are placeholder implementations for CI builds
// The actual implementations should be provided by the user

double ran1(long *idum) {
    // Simple linear congruential generator for testing
    static long seed = 12345;
    if (idum) seed = *idum;
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    return (double)seed / 2147483648.0;
}

double gasdev(long *idum) {
    // Simple Box-Muller transform for normal distribution
    static bool has_spare = false;
    static double spare;
    
    if (has_spare) {
        has_spare = false;
        return spare;
    }
    
    has_spare = true;
    double u = ran1(idum);
    double v = ran1(idum);
    double mag = sqrt(-2.0 * log(u));
    spare = mag * sin(2.0 * M_PI * v);
    return mag * cos(2.0 * M_PI * v);
}

double ran0(long *idum) {
    return ran1(idum);
}

double gammln(double xx) {
    // Simple approximation for log(gamma(x))
    return log(tgamma(xx));
}

double gammp(double a, double x) {
    // Incomplete gamma function P(a,x)
    // Simple approximation - not numerically accurate
    return 1.0 - exp(-x);
}

double gammq(double a, double x) {
    // Incomplete gamma function Q(a,x) = 1 - P(a,x)
    return exp(-x);
}