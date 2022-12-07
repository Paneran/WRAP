#include <complex.h>
#include <stdio.h>
#include "tools.h"
#define M_PI 3.1415926535897932384

int main() {
    double complex c2 = -1.0+4.0*I;
    double complex c1 = cexp(I*2*M_PI/2);
    double complex c3 = c1 * c2;
    
    double complex a[4] = {0, 1, 2, 3};
    double complex b[2] = {0, 1};

    double complex * A = FFT_r(a, 4);
    double complex * B = FFT_r(b, 2); 
    double complex * c = IFFT_r(A, 4);
    double complex * d = IFFT_r(B, 2);

    printf("Subtraction of complex numbers:\n");
    printf("values of complex number: A0=%.2lf+%.2lfi\n", creal(*(A)), cimag(*(A)));
    printf("values of complex number: A1=%.2lf+%.2lfi\n", creal(*(A+1)), cimag(*(A+1)));
    printf("values of complex number: A2=%.2lf+%.2lfi\n", creal(*(A+2)), cimag(*(A+2)));
    printf("values of complex number: A3=%.2lf+%.2lfi\n", creal(*(A+3)), cimag(*(A+3)));
    printf("values of complex number: B0=%.2lf+%.2lfi\n", creal(*(B)), cimag(*(B)));
    printf("values of complex number: B1=%.2lf+%.2lfi\n", creal(*(B+1)), cimag(*(B+1)));

    printf("values of complex number: c0=%.2lf+%.2lfi\n", creal(*(c)), cimag(*(c)));
    printf("values of complex number: c1=%.2lf+%.2lfi\n", creal(*(c+1)), cimag(*(c+1)));
    printf("values of complex number: c2=%.2lf+%.2lfi\n", creal(*(c+2)), cimag(*(c+2)));
    printf("values of complex number: c3=%.2lf+%.2lfi\n", creal(*(c+3)), cimag(*(c+3)));
    printf("values of complex number: d0=%.2lf+%.2lfi\n", creal(*(d)), cimag(*(d)));
    printf("values of complex number: d1=%.2lf+%.2lfi\n", creal(*(d+1)), cimag(*(d+1)));
}