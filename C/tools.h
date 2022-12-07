#include <complex.h>

// Main Functions
double complex* FFT_r(double complex* x, int N);
double complex* FFTdp(double* x, int N);
double * xcorr(double* x, int N);
double complex* IFFT_r(double complex* X, int N);
double* IFFTdp(double* X, int N);

// Helper Functions
int log2_(int N);