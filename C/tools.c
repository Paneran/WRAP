//==========================================
// Title:  Toolbox for WRAP MCU
// Author: Nathan Pereira
// Last Modified: December 7, 2022
//==========================================
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.1415926535897932384

const int max_N = 2048;
const int MAX_BUFFER = max_N*3 -2;

// log2 copied from Stack Overflow, by Leos313
int log2_(int N) {
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

// Radix-2 FFT with recursion. Solution Works, but is quite space inefficient.
// Requires N/2 recursive frames, creating arrays of length 1024 each instance.
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

// Radix-2 FFT with dynamic programming. Space complexity N.
// would it be worth it to make space NlogN to cut computation time in half?
void FFTdp(double complex* x, double complex* X, int N) {
    // allocate memory
    double complex X_e[max_N>>1];
    double complex X_o[max_N>>1];
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        *(X + i) = *(x + i);
    }
    
    // decompose x into its 2 point pairs for DFT. O(N log(N))
    int p = log2_(N);
    // Outside loop O(log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(X_e + k) = *(X + 2*k + j*(N>>i));
                *(X_o + k) = *(X + 2*k + 1 + j*(N>>i));
            }
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(X + k + j*(N>>i)) = *(X_e + k);
                *(X + k + j*(N>>i) + (N>>(i+1))) = *(X_o + k);
            }    
        }
    }
    
    // put together the 2 point FFT pairs to make whole FFT. 
    double complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = N>>(p-i); // half length of block
        w = cexp(-I*2*M_PI/(l<<1));
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                *(X_e + k) = *(X + j*(1<<(i+1)) + k);
                *(X_o + k) = *(X + j*(1<<(i+1)) + k + l);
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                *(X + k + j*(1<<(i+1))) = *(X_e + k) + *(X_o + k) * tw;
                *(X + k + j*(1<<(i+1)) + l) = *(X_e + k) - *(X_o + k) * tw;
                tw *= w;
            }
        }
    } 
}

// FFTdp and IFFTdp are nearly identical and easy to make one func
void IFFTdp(double complex* x, double complex* X, int N) {
    // allocate memory
    double complex X_e[max_N>>1];
    double complex X_o[max_N>>1];
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        *(X + i) = *(x + i);
    }
    
    // decompose x into its 2 point pairs for DFT. O(N log(N))
    int p = log2_(N);
    // Outside loop O(log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(X_e + k) = *(X + 2*k + j*(N>>i));
                *(X_o + k) = *(X + 2*k + 1 + j*(N>>i));
            }
            for (int k = 0; k < (N>>(i+1)); k++) {
                *(X + k + j*(N>>i)) = *(X_e + k);
                *(X + k + j*(N>>i) + (N>>(i+1))) = *(X_o + k);
            }    
        }
    }
    
    // put together the 2 point FFT pairs to make whole FFT. 
    double complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = N>>(p-i); // half length of block
        w = cexp(I*2*M_PI/(l<<1));
        double complex tw = 1;
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                *(X_e + k) = *(X + j*(1<<(i+1)) + k);
                *(X_o + k) = *(X + j*(1<<(i+1)) + l + k);
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                *(X + k + j*(1<<(i+1))) = *(X_e + k) + *(X_o + k) * tw;
                *(X + k + j*(1<<(i+1)) + l) = *(X_e + k) - *(X_o + k) * tw;
                tw *= w;
            }
        }
    } 

    // divide each element by N
    for (int i = 0; i < N; i++) {
        *(X+i) = *(X+i)/N;
    }
}

// crosscorrelation function
void xcorr(double* x1, double* x2, double* X, int l1, int l2) {
    // come up with a solution with pointers so we don't need to create a zero padded array
    // we can do better
    double x3[MAX_BUFFER];
    int z = 0, L;
    double * a, * b;
    if (l1 < l2) {
        z = l1 - 1;
        a = x2;
        b = x1;
        L = l2 + z<<1;
    } else {
        z = l2 - 1;
        a = x1;
        b = x2;
        L = l1 + z<<1;
    }

    for (int i = 0; i < L; i++) {
        *(X + i) = 0;
        if (i > z && i < L + z) {
            *(x3 + i) = *(a+i);
        } else {
            *(x3 + i) = 0;
        }
    }

    z++;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < z; j++) {
            *(X + i) += *(x3 + i + j) * *(x2 + j);
        }
    }
}

void complex2double(double complex* a, double * b) {
    return;
}

void double2complex(double complex* a, double * b) {
    return;
}
