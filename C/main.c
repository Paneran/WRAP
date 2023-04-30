#include <stdio.h>
#include "transmitter.h"
#include "receiver.h"
#include "constants.h"

int main() {
    // bits to send
    // packet header = {1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1};
    int key_w_bits[N] = {1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,
                         0,1,0,0,1,0,0,1,0,0,
                         1,0,0,0,0,0,0,1,1,0,
                         0,1,0,0,0,1,1,0,1,1,
                         1,1,0,1,1,0,1,1,1,0,
                         0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,1};
                   

    // initialize arrays and structs
    uint16_t output[(N)*(SPS)];
    params_t t_params = {.phase = 0};
    params_r r_params = {.CL_phase = 0, .TR_phase = 0, .sps = SPS};
    int symbs[N];
    uint8_t bits[N-45];
    int num_symbs = 0; 
    
    transmit(key_w_bits, output, &t_params);
    if(DEBUG) {
        FILE *fpt = fopen("data/MyFile.csv", "w+");
        int i = 0;
        for (; i < (N)*(SPS); i++) {
            fprintf(fpt, "%u,", output[i]);
        }
        fclose(fpt);
    }
    num_symbs = demodulate(output, symbs, &r_params);

    find_packet(symbs, bits, num_symbs);

    int sum = 0;

    for (int i = 0; i < N-45; i++) {
        printf("%i", bits[i]);
        sum += (bits[i]!=key_w_bits[i+45]);
    }
    printf("\n%i", sum);
    
}