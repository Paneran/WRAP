#include <stdint.h>
#include <math.h>
#include "tools.h"
#include "constants.h"
#include "transmitter.h"

// params
// output: output array
// bits: input bit array
// phase: phase to start modulator at 

void upsample(float * x, float * out, int size, int samppsymb);

void transmit(int * bits, uint16_t * output, params_t * params) {
    // create upsample array
    const int sample_sz = (N)*(SPS);
    float symbols[N];
    float samples[(N)*(SPS)];
    float out[1<<((int)ceil(log2((N)*(SPS))))];

    // create symbols
    for (int i = 0; i < N; i++) {
        symbols[i] = (bits[i]-0.5)*2;
    }
    
    // upsample bits
    upsample(symbols, samples, sample_sz, SPS);
    // convolve with pulse filter
    conv(samples, RRC, out, sample_sz, RRC_LEN);

    // shift and modulate
    float shift = RRC_LEN/2. - 0.5;
    int k;
    float dt = FC/FS;
    for (int i = shift ; i < sample_sz+RRC_LEN-1-shift; i++) {
        k = i - shift;
        out[k] = out[i]*cos(2*M_PI*dt*k + params->phase);
    }
    params->phase = 2*M_PI*dt*k + params->phase;

    // scale to output
    for (int i = 0; i < sample_sz; i++) {
        output[i] = (out[i]+1)*(1<<11);
    }
}

void upsample(float * x, float * out, int size, int samppsymb) {
    int k = 0;
    for (int i = 0; i < size; i++) {
        if (i%samppsymb == 0) {
            out[i] = x[k];
            k++;
        } else {
            out[i] = 0;
        }
    }
}