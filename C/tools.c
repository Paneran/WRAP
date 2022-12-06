#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.1415926535897932384

// Radix-2 FFT with recursion. Solution Works, but is quite space inefficient.
double complex* FFT_r(double* x, int N) {
    // Allocate memory for output FFT and split samples
    double complex* X = (double complex*) malloc(N*sizeof(double));
    double* x_e = (double*) malloc(N>>1 * sizeof(double));
    double* x_o = (double*) malloc(N>>1 * sizeof(double));
    // even and odd sampling
    for (int i = 0; i < (N>>1); i++) {
        *(x_e + i) = *(x + 2*i);
        *(x_o + i) = *(x + 2*i + 1);
    }
    // compute base case FFT
    if (N==2) {
        *X = *x_e + *x_o;
        *(X+1) = *x_e - *x_o;
        return X;
    }
    // compute FFT of substage
    double complex* X_e = FFT_r(x_e, N>>1);
    double complex* X_o = FFT_r(x_o, N>>1);
    // combine even and odd FFTs to make original FFT
    double complex w = cexp(-I*2*M_PI/N);
    double complex tw = 1;
    for (int k = 0; k < (N>>1); k++) {
        *(X + k) = *(X_e + k) + *(X_o + k) * tw;
        *(X + k + (N>>1)) = *(X_e + k) - *(X_o + k) * tw;
        tw = tw * w;
    }
    // return original FFT
    return X;
}

// Radix-2 FFT with dynamic programming. INCOMPLETE.
double complex* FFTdp(double* x, int N) {
    return NULL;
}