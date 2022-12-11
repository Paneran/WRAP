#include <complex.h>
#include <stdio.h>
#include "tools.h"

// testing program
int main() {
    int ft = 0;
    int cor = 0;
    int cast = 0;
    int fil = 1;
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
        double a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        double complex A[8];
        double b[2] = {0, 1};
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

    if (cast) {
        double a[3] = {3, 2, 7};
        double complex a_[3];
        double2complex(a, a_, 3);
        printf("a_:\n");
        for (int i = 0; i < 3; i++) {
            printf("values of complex number: %.2lf+%.2lfi\n", creal(*(a_+i)), cimag(*(a_+i)));
        }
        double b[3];
        complex2double(a_, b, 3);
        printf("b:\n");
        for (int i = 0; i < 3; i++) {
            printf("values of complex number: %.2lf\n", *(b+i));
        }
    }

    if (fil) {
        // fifth order low pass filter. generated with firpm
        // f = [0 0.2 0.25 1], a = [1, 1, 0, 0]: firpm(5, f, a)
        double lp[6] = {0.2630, 0.1926, 0.2207, 0.2207, 0.1926, 0.2630};
        double x[6] = {1, 1, 1, 1, 1, 1};
        double ans[6] = {0.2630, 0.4556, 0.6763, 0.8970, 1.0896, 1.3526};
        filter(x, lp, 6, 6);
        for (int i = 0; i < 6; i++) {
            printf("values of complex number: %.4lf\n", *(x+i));
        }
    }
}