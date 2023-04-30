#include <stdio.h>
#include "transmitter.h"
#include "receiver.h"
#include "constants.h"

int main() {
    // bits to send
    int bits[N] = {};

    // initialize arrays and structs
    uint16_t output[N*SPS];
    params_t t_params = {.phase = 0};
    params_r r_params = {.CL_phase = 0, .TR_phase = 0, .sps = SPS};
    int symbs[N];
    uint8_t out_bits[N];
    
    transmit(bits, output, &t_params);
    
    // print transmitted symbols to file
    FILE *fpt = fopen("MyFile.csv", "w+");
    int i = 0;
    for (; i < N*SPS; i++) {
        fprintf(fpt, "%u,", output[i]);
    }
    fclose(fpt);

    demodulate(output, symbs, &r_params);
    find_packet(symbs, bits, N);
    
}