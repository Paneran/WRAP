#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.1415926535897932384

int max_N = 2048;

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
void FFT_r(double complex* x, double complex * X, int N) {
    // Allocate memory for split samples
    double complex x_e[max_N>>1];
    double complex x_o[max_N>>1];
    // even and odd sampling
    for (int i = 0; i < (N>>1); i++) {
        *(x_e + i) = *(x + 2*i);
        *(x_o + i) = *(x + 2*i + 1);
    }
    // compute base case FFT
    if (N==2) {
        *X = *x_e + *x_o;
        *(X+1) = *x_e - *x_o;
        return;
    }
    // compute FFT of substage
    double complex X_e[max_N>>1]; 
    double complex X_o[max_N>>1]; 
    FFT_r(x_e, X_e, N>>1);
    FFT_r(x_o, X_o, N>>1);
    // combine even and odd sampled FFTs to make whole FFT
    double complex w = cexp(-I*2*M_PI/N);
    double complex tw = 1;
    for (int k = 0; k < (N>>1); k++) {
        *(X + k) = *(X_e + k) + *(X_o + k) * tw;
        *(X + (N>>1) + k) = *(X_e + k) - *(X_o + k) * tw;
        tw *= w;
    }
}

// Radix-2 FFT with dynamic programming. INCOMPLETE.
double complex* FFTdp(double complex* x, int N) {
    printf("b\n");
    // allocate memory
    double complex* X = (double complex*) malloc(N*sizeof(double));
    double complex* x_e = (double complex*) malloc(N>>1 * sizeof(double));
    double complex* x_o = (double complex*) malloc(N>>1 * sizeof(double));
    printf("b2\n");
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        *(X + i) = *(x + i);
    }
    
    // BELOW INCOMPLETE
    // break down x into smallest parts. O(N log(N))
    int p = log2_(N);
    // Outside loop O(log(N))
    printf("b3\n");
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        for (int j = 0; j < (1<<i); j++) {
            printf("%.2lf+%.2lfi\n", creal(*(X+1)), cimag(*(X+1)));
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(x_e + k) = *(X + 2*k + j*(N>>i));
                *(x_o + k) = *(X + 2*k + 1 + j*(N>>i));
            }
            printf("%.2lf+%.2lfi\n", creal(*(X)), cimag(*(X)));
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(X + k + j*(N>>i)) = *(x_e + k);
                *(X + k + j*(N>>i) + (N>>(i+1))) = *(x_o + k);
            }
            
        }
    }
    printf("you made it\n");
    return X;
}

void IFFT_r(double complex* X, double complex* x, int N) {
        // Allocate memory for output FFT and split samples
    double complex X_e[max_N>>1];
    double complex X_o[max_N>>1];
    // even and odd sampling
    for (int i = 0; i < (N>>1); i++) {
        *(X_e + i) = *(X + 2*i);
        *(X_o + i) = *(X + 2*i + 1);
    }
    // compute base case FFT
    if (N==2) {
        *x = (*X_e + *X_o)/2;
        *(x+1) = (*X_e - *X_o)/2;
        return;
    }
    // compute IFFT of substage
    double complex x_e[max_N>>1];
    double complex x_o[max_N>>1];
    IFFT_r(X_e, x_e, N>>1);
    IFFT_r(X_o, x_o, N>>1);
    // combine even and odd sampled FFTs to make whole FFT
    double complex w = cexp(I*2*M_PI/N);
    double complex tw = 1;
    for (int k = 0; k < (N>>1); k++) {
        *(x + k) = (*(x_e + k) + *(x_o + k) * tw)/2;
        *(x + (N>>1) + k) = (*(x_e + k) - *(x_o + k) * tw)/2;
        tw *= w;
    }
}

double* IFFTdp(double* X, int N) {
    return NULL;
}

// crosscorrelation function
double * xcorr(double* x, int N) {
    return NULL;
}
