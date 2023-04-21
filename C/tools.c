//==========================================
// Title:  Toolbox for WRAP MCU
// Author: Nathan Pereira
//==========================================
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ensure to set max_N to your max packet size. 
// keep this as a power of 2. 
#define max_N 32768

/* 
Radix-2 FFT with dynamic programming. Space complexity N.
Params:
x1: array to take FFT
X : output array
N: length of input array. Must be a power of 2.
*/
void FFT(const float * x, float complex* X, const int N) {
    int p = log2(N);
    if (floor(p) != p) {
        printf("WARNING: N must be a power of 2. Zero Pad first. it is %i\n", N);
        return; 
    }
    // allocate memory
    float complex X_e[max_N>>1];
    float complex X_o[max_N>>1];
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        X[i] = x[i];
    }
    
    // decompose x into its 2 point pairs for DFT. O(N log(N))
    // Outside loop O(log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        int b = N>>i;
        int l = b>>1;
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < l; k++) {
                X_e[k] = X[2*k + j*b];
                X_o[k] = X[2*k + 1 + j*b];
            }
            for (int k = 0; k < l; k++) {
                X[k + j*b] = X_e[k];
                X[k + j*b + l] = X_o[k];
            }    
        }
    }
    // put together the 2 point FFT pairs to make whole FFT. O(N log(N)) 
    float complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = 1<<i; // half block length
        int b = l<<1; // full block length
        w = cexp(-I*2*M_PI/(l<<1));
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                X_e[k] = X[j*b + k];
                X_o[k] = X[j*b + k + l];
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                X[k + j*b] = X_e[k] + X_o[k] * tw;
                X[k + j*b + l] = X_e[k] - X_o[k] * tw;
                tw *= w;
            }
        }
    } 
}

// FFT and IFFT are nearly identical and easy to make one func
/* 
Computes IFFT
Params:
x1: array to take IFFT
X : output array
N: length of input array. Must be a power of 2. 
*/
void IFFT(const float complex* x, float complex* X, const int N) {
    int p = log2(N);
    if (floor(p) != p) {
        printf("WARNING: N must be a power of 2. it is %i\n", N);
        return; 
    }
    // allocate memory
    float complex X_e[max_N>>1];
    float complex X_o[max_N>>1];
    // copy x into X so we don't modify x. 
    for (int i = 0; i < N; i++) {
        *(X + i) = *(x + i);
    }
    
    // decompose x into its 2 point pairs for IDFT. O(N log(N))
    // Outside loop O(log(N))
    for (int i = 0; i < p; i++) {
        // two for loops below O(N) time
        int b = N>>i;
        int l = b>>1;
        for (int j = 0; j < (1<<i); j++) {
            for (int k = 0; k < l; k++) {
                X_e[k] = X[2*k + j*b];
                X_o[k] = X[2*k + 1 + j*b];
            }
            for (int k = 0; k < l; k++) {
                X[k + j*b] = X_e[k];
                X[k + j*b + l] = X_o[k];
            }    
        }
    }
    // put together the 2 point IDFT pairs to make whole IFFT. O(N log(N)) 
    float complex w, tw;
    for (int i = 0; i < p; i++) {
        int l = 1<<i; // half block length
        int b = l<<1; // full block length
        w = cexp(I*2*M_PI/(l<<1));
        for (int j = 0; j < (N>>(i+1)); j++) { // iterate through each group
            for (int k = 0; k < l; k++) { // iterate through each point in each group
                // first half of group is X_e and second half is X_o
                X_e[k] = X[j*b + k];
                X_o[k] = X[j*b + k + l];
            }
            tw = 1;
            for (int k = 0; k < l; k++) {
                X[k + j*b] = X_e[k] + X_o[k] * tw;
                X[k + j*b + l] = X_e[k] - X_o[k] * tw;
                tw *= w;
            }
        }
    } 

    // divide each element by N
    for (int i = 0; i < N; i++) {
        X[i] = X[i]/N;
    }
}

/* 
crosscorrelation function
Params:
x1: array 1 to correlate
x2: array 2 to correlate
X : output array
l1: length of array 1
l2: length of array 2
L : length of array X
*/
void xcorr(float* x1, float* x2, float* X, const int l1, const int l2, const int L) {
    // preliminary warnings on output sizing. 
    if (L != (l1 + l2 -1)) {
        printf("WARNING: Make L=l1+l2-1. Op Cancelled.\n");
        return;
    }
    
    // set initial pointer positions of what values in each vector should be considered
    float * a1 = x1 + l1 - 1, *a2 = x1 + l1 - 1, *b1 = x2, *b2 = x2;
    // set all values in X to zero
    for (int i = 0; i < L; i++)
        *(X + i) = 0;
    
    // perform cross correlation
    int j;
    float * a, * b;
    for (int i = 1; i <= L; i++) {
        j = 0;
        // finding the correlation at a certain index. 
        while(1) {
            a = a2 + j;
            b = b1 + j;
            X[i-1] += *a * *b;
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

/* 
Convert complex array to float. Takes real component.
Params:
x1: array 1 to cast
X : output array
L : length of all arrays
*/
void complex2float(const float complex* x1, float * X, const int L) {
    for (int i = 0; i < L; i++) {
        X[i] = creal(x1[i]);
    }
}

/* 
Convert float array to complex
Params:
x1: array 1 to cast
X : output array
L : length of all arrays
*/
void float2complex(const float * x1, float complex * X, const int L) {
    for (int i = 0; i < L; i++) {
        X[i] = x1[i];
    }
}

/* 
Zeropads arrays to the nearest power of 2.
Pads >= l1 + l2 - 1 and to the nearest power of 2. 
if just want regular zero padding to power of 2, set l2 = 1
returns final size. 
Params:
x1: array to zeropad
X : output array
l1: length of input array
l2: length of another array(not passed in), which array should be padded to
*/
int zeropad(const float * x1, float * X, const int l1, const int l2) {
    int L = l1 + l2 - 1;
    int k = ceil(log2(L));
    int i = 0;
    for (; (i < 1<<k); i++) {
        if (i < l1) {
            X[i] = x1[i];
        } else {
            X[i] = 0;
        }
    }
    return i + 1;
}

/* 
Multiply complex float arrays together
Params:
x1: array 1 to multiply
x2: array 2 to multiply
X : output array
L : length of all arrays
*/
void multiply_c(const float complex *a, const float complex *b, float complex *c, int L) {
    for (int i = 0; i < L; i++) {
        c[i] = a[i] * b[i];
    }
}


/* 
Multiply float arrays together
Params:
x1: array 1 to multiply
x2: array 2 to multiply
X : output array
L : length of all arrays
*/
void multiply(const float *x1, const float *x2, float *X, int L) {
    for (int i = 0; i < L; i++) {
        X[i] = x1[i] * x2[i];
    }
}


/* 
filter takes "a" and modifies it with b
also must be of length 2^k where k is an integer
performs a linear convolution and stores the first l1 results back in a
have MATLAB generate the filters. 
Params:
x1: array 1 to filter
x2: array 2 to filter
l1: length of array 1
l2: length of array 2
output is stored in x1
*/
void filter(float * x1, const float * x2, const int l1, const int l2) {
    // could set a better upperbound/ design these functions s.t. no need for temp arrays
    float complex A[max_N];
    float complex B[max_N];
    float complex C[max_N];
    float temp[max_N];
    int L = zeropad(x1, temp, l1, l2);
    FFT(temp, A, L);
    zeropad(x2, temp, l1, l2);
    FFT(temp, B, L);
    multiply_c(A, B, C, L);
    IFFT(C, A, L);
    complex2float(A, x1, l1);
}


/* 
regular linear convolution
Params:
x1: array 1 to convolve
x2: array 2 to convolve
X : output array
l1: length of array 1
l2: length of array 2
L : length of array X
*/
void conv(const float * x1, const float * x2, float * X, const int l1, const int l2, const int L) {
    if (L != (l1 + l2 -1)) {
        printf("WARNING: Make L=l1+l2-1. Op Cancelled.\n");
        return;
    }
    int N = 1<<((int) ceil(log2(L)));
    float complex A[max_N];
    float complex B[max_N];
    float complex C[max_N];
    float temp[max_N];
    zeropad(x1, temp, l1, l2);
    // failing at FFT for some reason
    FFT(temp, A, N);
    zeropad(x2, temp, l2, l1);
    FFT(temp, B, N);
    multiply_c(A, B, C, N);
    IFFT(C, A, N);
    complex2float(A, X, L);
} 