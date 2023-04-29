#include <stdio.h>
#include "transmitter.h"
#include "receiver.h"
#include "constants.h"

int main() {
    // bits to send
    int key_w_bits[N] = {1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1,
                     1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1,
                     1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0,1,0,0,1,0,0,1,0,0,
                   1,0,0,0,0,0,0,1,1,0,
                   0,1,0,0,0,1,1,0,1,1,
                   1,1,0,1,1,0,1,1,1,0,
                   0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,
        0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,
        1,1,0,1,0,0,0,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,
        0,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,
        0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,1};

    // initialize arrays and structs
    uint16_t output[N*SPS];
    params_t t_params = {.phase = 0};
    params_r r_params = {.CL_phase = 0, .TR_phase = 0, .sps = SPS};
    int symbs[N];
    uint8_t bits[N];
    
    transmit(key_w_bits, output, &t_params);

    demodulate(output, symbs, &r_params);
    find_packet(symbs, bits, N);


    
    FILE *fpt = fopen("MyFile.csv", "w+");
    int i = 0;
    for (; i < N; i++) {
        fprintf(fpt, "%u,", bits[i]);
    }
    fclose(fpt);
    
}