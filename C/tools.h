#ifndef TOOLS
#define TOOLS

#include <complex.h>
#include <stdlib.h>

// Main Functions
void xcorr(double* x1, double* x2, double* X, const int l1, const int l2, const int L);
void conv(const double * x1, const double * x2, double * X, const int l1, const int l2, const int L);
void filter(double * x1, const double * x2, const int l1, const int l2);

void complex2double(const double complex* x1, double * X, const int L);
void double2complex(const double * x1, double complex * X, const int L);
void multiply_c(const double complex *a, const double complex *b, double complex *c, int L);
void multiply(const double *x1, const double *x2, double *X, int L);
void FFT(const double * x, double complex* X, const int N);
void IFFT(const double complex* x, double complex* X, const int N);
int zeropad(const double * x1, double * X, const int l1, const int l2);

#endif