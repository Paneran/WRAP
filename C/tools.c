#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.1415926535897932384

int log2_(int N) {
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

// Radix-2 FFT with recursion. Solution Works, but is quite space inefficient.
// Requires N/2 recursive frames...
// need to free up memory
double complex* FFT_r(double complex* x, int N) {
    // Allocate memory for output FFT and split samples
    double complex* X = (double complex*) malloc(N*sizeof(double));
    double complex* x_e = (double complex*) malloc(N>>1 * sizeof(double));
    double complex* x_o = (double complex*) malloc(N>>1 * sizeof(double));
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
    // combine even and odd sampled FFTs to make whole FFT
    double complex w = cexp(-I*2*M_PI/N);
    double complex tw = 1;
    for (int k = 0; k < (N>>1); k++) {
        *(X + k) = *(X_e + k) + *(X_o + k) * tw;
        *(X + (N>>1) + k) = *(X_e + k) - *(X_o + k) * tw;
        tw *= w;
    }
    // return FFT
    return X;
}

// Radix-2 FFT with dynamic programming. INCOMPLETE.
double complex* FFTdp(double* x, int N) {
    // allocate memory
    double complex* X = (double complex*) malloc(N*sizeof(double));
    double* x_e = (double*) malloc(N>>1 * sizeof(double));
    double* x_o = (double*) malloc(N>>1 * sizeof(double));
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        *(X + i) = *(x + i);
    }

    // BELOW INCOMPLETE
    // break down x into smallest parts. O(N log(N))
    int p = log2_(N);
    // Outside loop O(log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < (N>>i); k++) {
                *(x_e + k) = *(X + 2*k + j*(N>>i));
                *(x_o + k) = *(X + 2*k + 1 + j*(N>>i));
            }
            for (int k = 0; k < (N>>i); k++) {
                *(X + k + j*(N>>i)) = *(x_e + k);
                *(X + k + j*(N>>i) + (N>>(i+1))) = *(x_o + k);
            }
        }
    }
    return X;
}

double complex* IFFT_r(double complex* X, int N) {
        // Allocate memory for output FFT and split samples
    double complex* x = (double complex*) malloc(N*sizeof(double));
    double complex* X_e = (double complex*) malloc(N>>1 * sizeof(double));
    double complex* X_o = (double complex*) malloc(N>>1 * sizeof(double));
    // even and odd sampling
    for (int i = 0; i < (N>>1); i++) {
        *(X_e + i) = *(X + 2*i);
        *(X_o + i) = *(X + 2*i + 1);
    }
    // compute base case FFT
    if (N==2) {
        *x = *X_e + *X_o;
        *(x+1) = *X_e - *X_o;
        return X;
    }
    // compute IFFT of substage
    double complex* x_e = IFFT_r(X_e, N>>1);
    double complex* x_o = IFFT_r(X_o, N>>1);
    // combine even and odd sampled FFTs to make whole FFT
    double complex w = cexp(I*2*M_PI/N);
    double complex tw = 1;
    for (int k = 0; k < (N>>1); k++) {
        *(x + k) = *(x_e + k) + *(x_o + k) * tw;
        *(x + (N>>1) + k) = *(x_e + k) - *(x_o + k) * tw;
        tw *= w;
    }
    // return FFT
    return x;
}

double* IFFTdp(double* X, int N) {
    return NULL;
}

// crosscorrelation function
double * xcorr(double* x, int N) {
    return NULL;
}
