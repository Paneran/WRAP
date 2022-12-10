#include <complex.h>
#include <stdlib.h>

// Main Functions
void FFTdp(double * x, double complex * X, int N);
void IFFTdp(double complex * X, double complex * x, int N);
void xcorr(double* x1, double* x2, double* X, int l1, int l2, int L);
void complex2double(double complex* a, double * b, int L);
void double2complex(double * a, double complex * b, int L);
void filter(double * a, double * b, int l1, int l2);
void multiply_c(complex double *a, complex double *b, complex double *c, int L);
void multiply(double *a, double *b, double *c, int L);
int zeropad(double * a, double * b, int l1, int l2);

// Helper Functions
int log2_(int N);
