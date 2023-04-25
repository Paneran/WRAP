#ifndef RECEIVE
#define RECEIVE

int demodulate(uint16_t * samples, int * symbs, struct parameters_r * params);
void find_packet(int * symbs, uint8_t * out, int num_symbs);

#endif