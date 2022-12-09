//==========================================
// Title:  Toolbox for WRAP MCU
// Author: Nathan Pereira
// Last Modified: December 8, 2022
//==========================================
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.1415926535897932384

const int max_N = 2048;

// log2 copied from Stack Overflow, by Leos313
int log2_(int N) {
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
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
        int b = N>>i;
        int l = b>>1;
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < l; k++) {
                *(X_e + k) = *(X + 2*k + j*b);
                *(X_o + k) = *(X + 2*k + 1 + j*b);
            }
            for (int k = 0; k < l; k++) {
                *(X + k + j*b) = *(X_e + k);
                *(X + k + j*b + l) = *(X_o + k);
            }    
        }
    }
    
    // put together the 2 point FFT pairs to make whole FFT. O(N log(N)) 
    double complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = 1<<i; // half block length
        int b = l<<1; // full block length
        w = cexp(-I*2*M_PI/(l<<1));
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                *(X_e + k) = *(X + j*b + k);
                *(X_o + k) = *(X + j*b + k + l);
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                *(X + k + j*b) = *(X_e + k) + *(X_o + k) * tw;
                *(X + k + j*b + l) = *(X_e + k) - *(X_o + k) * tw;
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
    
    int p = log2_(N);
    // decompose x into its 2 point pairs for IDFT. O(N log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        int b = N>>i;
        int l = b>>1;
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < l; k++) {
                *(X_e + k) = *(X + 2*k + j*b);
                *(X_o + k) = *(X + 2*k + 1 + j*b);
            }
            for (int k = 0; k < l; k++) {
                *(X + k + j*b) = *(X_e + k);
                *(X + k + j*b + l) = *(X_o + k);
            }    
        }
    }
    
    // put together the 2 point IDFT pairs to make whole IFFT. O(N log(N)) 
    double complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = 1<<i; // half block length
        int b = l<<1; // full block length
        w = cexp(I*2*M_PI/(l<<1));
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                *(X_e + k) = *(X + j*b + k);
                *(X_o + k) = *(X + j*b + k + l);
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                *(X + k + j*b) = *(X_e + k) + *(X_o + k) * tw;
                *(X + k + j*b + l) = *(X_e + k) - *(X_o + k) * tw;
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
void xcorr(double* x1, double* x2, double* X, int l1, int l2, int L) {
    // preliminary warnings on output sizing. 
    if (L < (l1 + l2 -1)) {
        printf("WARNING: seg fault may occur. Make L=l1+l2-1\n");
        L = l1 + l2 - 1;
    } else if (L > (l1 + l2 - 1)) {
        printf("WARNING: array size bigger than it needs to be. Make L=l1+l2-1\n");
        L = l1 + l2 - 1;
    }
    
    // set initial pointer positions of what values in each vector should be considered
    double * a1 = x1 + l1 - 1, *a2 = x1 + l1 - 1, *b1 = x2, *b2 = x2;
    // set all values in X to zero
    for (int i = 0; i < L; i++)
        *(X + i) = 0;
    
    // perform cross correlation
    int j;
    double * a, * b;
    for (int i = 1; i <= L; i++) {
        j = 0;
        // finding the correlation at a certain index. 
        while(1) {
            a = a2 + j;
            b = b1 + j;
            *(X + i-1) += *a * *b;
            if (a == a1 || b == b2)
                break;
            j++;
        }
        // conditional statements to move boundary pointers. 
        (i>= l2) ? (a1--) : (0);
        (i < l1) ? (a2--) : (0);
        (i>= l1) ? (b1++) : (0);
        (i < l2) ? (b2++) : (0);
    }
}

void complex2double(double complex* a, double * b) {
    return;
}

void double2complex(double complex* a, double * b) {
    return;
}

// filter takes "a" and modifies it with b. 
// have MATLAB generate the filters. 
void filter(double * a, double * b, int l1, int l2) {
    return;
}
