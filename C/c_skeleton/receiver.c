#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "tools.h"
#include "constants.h"
#include "receiver.h"

void normalize(const uint16_t * samples, float * norm_samples);
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

int demodulate(const uint16_t * samples, int * symbs, params_r * params) {
    // normalize the input
    
    // Costas Loop
    
    // filter w SRRC

    // calculate zero crossings
    
    // timing recovery
    int bit_len = 0; // number of samples recovered from timing recovery

    return bit_len;
}


void find_packet(int * symbs, uint8_t * out, const int num_symbs) {
    // take cross correlation
    
    // find packet
    
    // convert symbols to bits

}