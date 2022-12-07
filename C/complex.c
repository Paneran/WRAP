#include <complex.h>
#include <stdio.h>
#include "tools.h"
#define M_PI 3.1415926535897932384

int main() {
    double complex c2 = -1.0+4.0*I;
    double complex c1 = cexp(I*2*M_PI/2);
    double complex c3 = c1 * c2;

    double complex a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    double complex A[8];
    double complex b[2] = {0, 1};
    double complex B[2];
    double complex c[8]; 
    double complex d[8];

    FFT_r(a, A, 8);
    FFT_r(b, B, 2); 
    IFFT_r(A, c, 8);
    IFFT_r(B, d, 2);

    
    printf("A:\n");
    for (int i = 0; i < 8; i++) {
        printf("values of complex number: %.2lf+%.2lfi\n", creal(*(A+i)), cimag(*(A+i)));
    }
    printf("B:\n");
    for (int i = 0; i < 2; i++) {
        printf("values of complex number: %.2lf+%.2lfi\n", creal(*(B+i)), cimag(*(B+i)));
    }
    printf("c:\n");
    for (int i = 0; i < 8; i++) {
        printf("values of complex number: %.2lf+%.2lfi\n", creal(*(c+i)), cimag(*(c+i)));
    }
    printf("d:\n");
    for (int i = 0; i < 2; i++) {
        printf("values of complex number: %.2lf+%.2lfi\n", creal(*(d+i)), cimag(*(d+i)));
    }
}