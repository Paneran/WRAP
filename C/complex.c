#include <complex.h>
#include <stdio.h>
#define M_PI 3.1415926535897932384

int main() {
    double complex c2 = -1.0+4.0*I;
    double complex c1 = cexp(I*2*M_PI/2);
    double complex c3 = c1 * c2;
    
    printf("Subtraction of complex numbers:\n");
    printf("values of complex number c1: c1=%.2lf+%.2lfi\n", creal(c1), cimag(c1));
}