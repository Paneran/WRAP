#ifndef TOOLS
#define TOOLS

#include <complex.h>
#include <stdlib.h>

// Main Functions
void xcorr(int * x1, int * x2, int * X, const int l1, const int l2);
void conv(const float * x1, const float * x2, float * X, const int l1, const int l2);
void filter(float * x1, const float * x2, float * out, const int l1, const int l2);

void complex2float(const float complex* x1, float * X, const int L);
void float2complex(const float * x1, float complex * X, const int L);
void multiply_c(const float complex *a, const float complex *b, float complex *c, int L);
void multiply(const float *x1, const float *x2, float *X, int L);
void FFT(const float * x, float complex* X, const int N);
void IFFT(const float complex* x, float complex* X, const int N);
int zeropad(const float * x1, float * X, const int l1, const int l2);
float wrap_to_pi(const float x);

#endif