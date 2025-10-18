#ifndef ZROOTS2_H
#define ZROOTS2_H

#include <complex>

// Stub header for numerical recipes zroots functions
// These are placeholder declarations for CI builds
// The actual implementations should be provided by the user

void zroots(double a[], int m, double roots[], bool polish, const char* name = nullptr);
void zroots(std::complex<double> a[], int m, std::complex<double> roots[], bool polish, const char* name = nullptr);

#endif // ZROOTS2_H
