#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "tools.h"
#include "constants.h"

const float lp[ORDER + 1] = {};
int key[15] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1}; // packet header
const float RRC_r[RRC_LEN] = {};

void costas_loop(float * samples, float * samples_d);
int  timing_recovery(float * filtered_samps, int * symbs);

/*
Takes samples and turns them into symbols.
samples: samples from ADC
symbs: symbs at the output
params: parameters that need to be stored intermediately
returns length of symbs array. symbs array must be 
allocated for longer than samples/sps + some margin
*/

int demodulate(uint16_t * samples, int * symbs, params_r * params) {
    // normalize samples
    
    // costas loop

    // filter w SRRC
    
    // calculate zero crossings

    // timing recovery

    return 0;
}


void find_packet(int * symbs, uint8_t * out, int num_symbs) {
    // take cross correlation
    
    // find packet
    
    // convert symbols to bits
    
}