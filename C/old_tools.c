#include <complex.h>

// Radix-2 FFT with recursion. Solution Works, but is quite space inefficient.
// Requires N/2 recursive frames, creating arrays of length 1024 each instance.
// USE THE DYNAMIC PROGRAMMING FFTs/IFFTs
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
