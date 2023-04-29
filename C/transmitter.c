#include <stdint.h>
#include <math.h>
#include "tools.h"
#include "constants.h"

// params
// output: output array
// bits: input bit array
// phase: phase to start modulator at 

float F = 1; // Fc/Fs
const float RRC_t[RRC_LEN] = {0.0290581511866829,0.0292861195166709,0.0287874417058172,0.0275436194555175,0.0255527588425467,0.0228302397253860,0.0194090371891105,0.0153396797427319,
                             0.0106898363592876,0.00554353211700836,-8.89649295997777e-18,-0.00582781581531654,-0.0118150822918442,-0.0178271953766885,-0.0237221565644684,
                             -0.0293531653612106,-0.0345713796105044,-0.0392287913457370,-0.0431811625587258,-0.0462909631069959,-0.0484302519778048,-0.0494834433212715,
                             -0.0493499000671152,-0.0479463005336786,-0.0452087271829673,-0.0410944315056948,-0.0355832348467026,-0.0286785316929336,-0.0204078694131854,
                             -0.0108230865141592,8.89649295997777e-18,0.0119623587788076,0.0249429515050043,0.0388003664080867,0.0533748522700539,0.0684907191761580,
                             0.0839590647683679,0.0995807780314863,0.115149766823269,0.130456350574261,0.145290755933415,0.159446650701875,0.172724650234903,0.184935730629903,
                             0.195904484459525,0.205472157528474,0.213499409080215,0.219868742979158,0.224486563545039,0.227284816797345,0.228222185737461,0.227284816797345,
                             0.224486563545039,0.219868742979158,0.213499409080215,0.205472157528474,0.195904484459525,0.184935730629903,0.172724650234903,0.159446650701875,
                             0.145290755933415,0.130456350574261,0.115149766823269,0.0995807780314863,0.0839590647683679,0.0684907191761580,0.0533748522700539,0.0388003664080867,
                             0.0249429515050043,0.0119623587788076,8.89649295997777e-18,-0.0108230865141592,-0.0204078694131854,-0.0286785316929336,-0.0355832348467026,
                             -0.0410944315056948,-0.0452087271829673,-0.0479463005336786,-0.0493499000671152,-0.0494834433212715,-0.0484302519778048,-0.0462909631069959,
                             -0.0431811625587258,-0.0392287913457370,-0.0345713796105044,-0.0293531653612106,-0.0237221565644684,-0.0178271953766885,-0.0118150822918442,
                             -0.00582781581531654,-8.89649295997777e-18,0.00554353211700836,0.0106898363592876,0.0153396797427319,0.0194090371891105,0.0228302397253860,
                             0.0255527588425467,0.0275436194555175,0.0287874417058172,0.0292861195166709,0.0290581511866829};

void upsample(float * x, float * out, int size, int samppsymb);

void transmit(int * bits, uint16_t * output, params_t * params) {
    // create upsample array
    const int sample_sz = N*SPS;
    float symbols[N];
    float samples[N*SPS];
    float out[1<<((int)ceil(log2(N*SPS)))];
    float t[N*SPS];

    // create symbols
    for (int i = 0; i < N; i++) {
        symbols[i] = (bits[i]-0.5)*2;
    }
    
    // upsample bits
    upsample(symbols, samples, sample_sz, SPS);
    // convolve with pulse filter
    conv(samples, RRC_t, out, sample_sz, RRC_LEN);
    // narrow the range of convolution while also finding min and max
    float shift = RRC_LEN/2. - 0.5;
    int k;
    for (int i = shift ; i < sample_sz+RRC_LEN-1-shift; i++) {
        k = i - shift;
        out[k] = out[i]*cos(2*M_PI*k*F + params->phase);
    }
    params->phase = 2*M_PI*k*F + params->phase;

    // scale to output
    for (int i = 0; i < sample_sz; i++) {
        output[i] = (out[i]+1)*(1<<11);
    }
}

void upsample(float * x, float * out, int size, int samppsymb) {
    int k = 0;
    for (int i = 0; i < size; i++) {
        if (i%samppsymb == 0) {
            out[i] = x[k];
            k++;
        } else {
            out[i] = 0;
        }
    }
}