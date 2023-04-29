#ifndef RECEIVE
#define RECEIVE

#include "constants.h"

int demodulate(uint16_t * samples, int * symbs, params_r * params);
void find_packet(int * symbs, uint8_t * out, int num_symbs);

#endif