#include <complex.h>
#include <stdlib.h>

// Main Functions
void xcorr(double* x1, double* x2, double* X, int l1, int l2, int L);
void conv(double * x1, double * x2, double * X, int l1, int l2, int L);
void filter(double * a, double * b, int l1, int l2);

void complex2double(double complex * a, double * b, int L);
void double2complex(double * a, double complex * b, int L);
void multiply_c(double complex *a, double complex * b, double complex *c, int L);
void multiply(double *a, double * b, double *c, int L);
int zeropad(double * a, double * b, int l1, int l2);
void FFT(double * x, double complex * X, int N);
void IFFT(double complex * X, double complex * x, int N);
