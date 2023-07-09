#ifndef RECEIVE
#define RECEIVE

#include "constants.h"

#define ORDER 5     // low pass filter order
#define CORRELATION_BUFFER 32768 // length of how long the correlation is
#define NUM_SAMPLES (N*SPS)       // Number of samples to process in receiver buffer (length of buffer) 

// parameters to pass between calls of receiver
typedef struct params_r {
    float CL_phase; // costas loop phase
    float CL_integrator;
    float TR_phase; // timing recovery phase
    float TR_integrator;
    float sps;      // last sps
} params_r;

static const float lp[ORDER + 1] = {0.238053579626575,0.064423148855147,0.0684583456734438,0.0684583456734438,0.0644231488551473,0.238053579626575};
static int key[15] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1}; // packet header

int demodulate(const uint16_t * samples, int * symbs, params_r * params);
uint8_t find_packet(int * symbs, uint8_t * out, const int num_symbs);

#endif