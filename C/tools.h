#include <complex.h>

// Main Functions
void FFT_r(double complex* x, double complex * X, int N);
double complex* FFTdp(double complex* x, int N);
double * xcorr(double* x, int N);
void IFFT_r(double complex* X, double complex* x, int N);
double* IFFTdp(double* X, int N);

// Helper Functions
int log2_(int N);