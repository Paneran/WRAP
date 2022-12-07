#include <complex.h>

// Main Functions
void FFT_r(double complex* x, double complex* X, int N);
void FFTdp(double complex* x, double complex* X, int N);
void IFFT_r(double complex* X, double complex* x, int N);
void IFFTdp(double complex* X, double complex* x, int N);
void xcorr(double* x, double* X, int N);

// Helper Functions
int log2_(int N);