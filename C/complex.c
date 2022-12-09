#include <complex.h>
#include <stdio.h>
#include "tools.h"
#define M_PI 3.1415926535897932384

int main() {
    int ft = 1;
    int cor = 0;
    if (cor) {
        double e[8] = {1, 2, 3, 4, 6, 9, 2, 4};
        double f[8] = {2, 3, 5, 7, 0, 2, 9, 5};
        double g[15];
        xcorr(f, e, g, 8, 8, 15);
        printf("g:\n");
        for (int i = 0; i < 15; i++) {
            printf("values of number: %.2lf\n", *(g+i));
        }

        xcorr(e, f, g, 8, 8, 15);
        printf("g:\n");
        for (int i = 0; i < 15; i++) {
            printf("values of number: %.2lf\n", *(g+i));
        }
    }
    

    if (ft) {
        double complex a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        double complex A[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        double complex b[2] = {0, 1};
        double complex B[2];
        double complex c[8]; 
        double complex d[8];

        FFTdp(a, A, 8);
        FFTdp(b, B, 2); 
        IFFTdp(A, c, 8);
        IFFTdp(B, d, 2);

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
}