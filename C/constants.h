#ifndef CONSTANTS
#define CONSTANTS

// ***** CHANGE ALL PARAMS ACCORDING TO WHAT YOU USE ON YOUR CODE *****

#define RRC_LEN 101             // length of square root raised cosine filter
#define N 301                   // number of symbols being transmitted (data and header)
#define SPS 100                 // sps

// receiver params
#define ORDER 5                 // low pass filter order
#define CORRELATION_BUFFER 1024 // length of how long the correlation is
#define NUM_SAMPLES 30100       // Number of samples to process in receiver buffer (length of buffer) 

// parameters to pass between calls of transmitter
typedef struct params_t {
    float phase; // initial phase of modulating signal
} params_t;

// parameters to pass between calls of receiver
typedef struct params_r {
    float CL_phase; // costas loop phase
    float TR_phase; // timing recovery phase
    float sps;      // last sps
} params_r;

#endif