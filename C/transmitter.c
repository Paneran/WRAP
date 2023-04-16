#include <complex.h>
#include "tools.h"
#include "transmitter.h"


// params
// output: output array
// bits: input bit array
// N: length of bits array
#define N 256
#define sps 100

void calculate_transmit(uint16_t * output, double * bits) {
    // create upsample array
    const int sample_sz = (N+45)*sps;
    double complex samples[(N+45)*sps];
    
    // upsample bits
    int k = 0;
    for (int i = 0; i < sample_sz; i++) {
        if (i%sps == 0) {
            if (k < 45) {
                samples[i] = key[k];
            } else {
                samples[i] = bits[k];
            }
            k++;
        } else {
            samples[i] = 0;
        }
    }

    

    




    
    
}