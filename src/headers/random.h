#ifndef RANDOM_H
#define RANDOM_H

// Stub header for numerical recipes random functions
// These are placeholder declarations for CI builds
// The actual implementations should be provided by the user

double ran1(long *idum);
double ran2(long *idum);
double gasdev(long *idum);
double ran0(long *idum);
double gammln(double xx);
double gammp(double a, double x);
double gammq(double a, double x);
int randint(int min, int max, long *seed);
double poisson(double mean, long *seed);

#endif // RANDOM_H
