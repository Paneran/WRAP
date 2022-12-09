#include <complex.h>
#include <stdlib.h>

// Main Functions
void FFTdp(double complex * x, double complex * X, int N);
void IFFTdp(double complex * X, double complex * x, int N);
void xcorr(double* x1, double* x2, double* X, int l1, int l2, int L);

// Helper Functions
int log2_(int N);