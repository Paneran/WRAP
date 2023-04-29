#include <stdint.h>
#include <math.h>
#include "tools.h"
#include "constants.h"

// params
// output: output array
// bits: input bit array
// phase: phase to start modulator at 

float F = 1; // Fc/Fs
const float RRC_t[RRC_LEN] = {};

void transmit(int * bits, uint16_t * output, params_t * params) {
    // create upsample array

    // create symbols

    // upsample symbols

    // convolve with pulse filter
    
    // scale to output

}
