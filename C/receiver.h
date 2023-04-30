#ifndef RECEIVE
#define RECEIVE

#include "constants.h"

int demodulate(const uint16_t * samples, int * symbs, params_r * params);
void find_packet(int * symbs, uint8_t * out, const int num_symbs);

#endif