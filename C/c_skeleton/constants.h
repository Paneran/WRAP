#ifndef CONSTANTS
#define CONSTANTS

// ***** CHANGE RRC AND N AS YOU SEE FIT *****

#define DEBUG 1         // if want to print .csv files

#define RRC_LEN 401     // length of square root raised cosine filter
#define N 256+45        // number of symbols being transmitted (data and header)
#define FS 4000000      // sampling frequency
#define FC 1000000.     // carrier frequency
#define RS 50000        // symbol rate
#define SPS FS/RS       // samples per symbol

static const float SRRC[RRC_LEN] ={}; // SRRC filter

#endif