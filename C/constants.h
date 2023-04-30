#ifndef CONSTANTS
#define CONSTANTS

// ***** CHANGE ALL PARAMS ACCORDING TO WHAT YOU USE ON YOUR CODE *****

#define DEBUG 1

#define RRC_LEN 401             // length of square root raised cosine filter
#define N 256                   // number of symbols being transmitted (data and header)
#define FS 4000000
#define FC 1000000.
#define RS 50000
#define SPS FS/RS 

// receiver params
#define ORDER 5                 // low pass filter order
#define CORRELATION_BUFFER 32768 // length of how long the correlation is
#define NUM_SAMPLES N*SPS       // Number of samples to process in receiver buffer (length of buffer) 

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