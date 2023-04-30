#ifndef TRANSMIT
#define TRANSMIT

#include <stdint.h>
#include "constants.h"

// parameters to pass between calls of transmitter
typedef struct params_t {
    float phase; // initial phase of modulating signal
} params_t;

void transmit(int * bits, uint16_t * output, params_t * params);

#endif