#include "tools.h"
#include <math.h>

#define ORDER 5
#define NUM_SAMPLES 30000 // ***** change *****
#define RRC_len 101
#define CORRELATION_BUFFER 1024
#define INFO_SIZE 256

struct parameters_r {
    float phase;
    float sps;
};
void normalize(uint16_t * samples, float * norm_samples);
void costas_loop(float * samples, float * samples_d);

const int key[15] = {1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1};
const float lp[ORDER + 1] = {0.238053579626575,0.064423148855147,0.0684583456734438,0.0684583456734438,0.0644231488551473,0.238053579626575};
const float RRC[RRC_len] = {0.0290581511866829,0.0292861195166709,0.0287874417058172,0.0275436194555175,0.0255527588425467,0.0228302397253860,0.0194090371891105,0.0153396797427319,
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

/*
Takes samples and turns them into symbols.
samples: samples from ADC
symbs: symbs at the output
params: parameters that need to be stored intermediately
returns length of symbs array. symbs array must be 
allocated for longer than samples/sps + some margin
*/

int demodulate(uint16_t * samples, int * symbs, struct parameters_r * params) {
    float norm_samples[NUM_SAMPLES];
    normalize(samples, norm_samples);

    float samples_d[NUM_SAMPLES] = {0, 0, 0, 0, 0, 0};
    costas_loop(norm_samples, samples_d);
    
    // filter w SRRC
    float filtered_samps[1<<((int)ceil(log2(NUM_SAMPLES)))];
    conv(samples_d, RRC, filtered_samps, NUM_SAMPLES, RRC_len);

    // timing recovery
    return timing_recovery(filtered_samps, symbs);
}


void find_packet(int * symbs, uint8_t * out, int num_symbs) {
    // take cross correlation
    float xcorr_out[CORRELATION_BUFFER];
    xcorr(symbs, key, xcorr_out, num_symbs, 15);

    // find packet
    int shift = 0;
    for (int i = 0; i < num_symbs + 14; i++) {
        if (abs(xcorr_out[i] + xcorr_out[i+15] + xcorr_out[i+30] ) > 39) {
            shift = i+45;
            if (xcorr_out[i] < 0) {
                for (int i = 0; i < INFO_SIZE; i++) {
                    symbs[shift + i] = symbs[shift+ i]*-1;
                }
            }
            break;
        }
    }
    
    // convert symbols to bits
    for (int i = 0; i < INFO_SIZE; i++) {
        out[i] = (symbs[shift+i]+1)*0.5;
    }
}

void normalize(uint16_t * samples, float * norm_samples) {
        // Normalize signal
    float var = 0, mean = 0;
    // find mean
    for (int i = 0; i < NUM_SAMPLES; i++) {
        mean += samples[i];
    }
    mean /= NUM_SAMPLES;
    // find sample variance
    for (int i = 0; i < NUM_SAMPLES; i++) {
        float temp = samples[i]-mean;
        var += temp * temp;
    }
    var /= NUM_SAMPLES-1;

    // normalize
    for (int i = 0; i < NUM_SAMPLES; i++) {
        norm_samples[i] = (samples[i]-mean)/var;
    }
}

void costas_loop(float * samples, float * samples_d) {
    // Costas Loop
    const int f_0 = 1000000;
    float t = 0;
    float phase = 0;
    float inph[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float quad[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float inph_[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float quad_[ORDER+1] = {0, 0, 0, 0, 0, 0};
    double error = 0;
    double integrator = 0;

    float kp = 0.05;
    float ki = 0.0005;

    for (int i = ORDER; i < NUM_SAMPLES+ORDER; i++) {
        // define t from microcontroller
        inph_[ORDER] = samples[i-ORDER]*2*cos(2*M_PI*f_0*t + phase);
        quad_[ORDER] = samples[i-ORDER]*-2*sin(2*M_PI*f_0*t + phase);

        filter(inph_,lp, inph, ORDER+1,ORDER+1);
        filter(quad_,lp, quad, ORDER+1,ORDER+1);

        samples_d[i] = inph[ORDER];

        error = inph[ORDER] * quad[ORDER];
        integrator = integrator + ki*error;
        phase = phase + kp*error + integrator;
        
        // shift the values of inph_ and quad_
        for (int j = 1; j < ORDER; j++) {
            inph_[j-1] = inph_[j];
            quad_[j-1] = quad_[j];
        }
    }
}

int timing_recovery(float * filtered_samps, int * symbs) {
    const float kp_PLL = 0.3;
    const float ki_PLL = 1.1;
    const int guess_sps = 20;
    const float margin = 0.75;

    float sps = params.sps;
    float integrator = 0, error = 0;
    float phase, prev_phase = params.phase;
    int bit_len = 0;

    float zc[NUM_SAMPLES-1];
    // calculate zero crossings
    for (int i = 0; i < NUM_SAMPLES-1; i++) {
        zc[i] = !(sign(filtered_samps[i+1] * filtered_samps[i])+1);
    }

    // timing recovery
    for (int i = 1; i < NUM_SAMPLES; i++) {
        phase = prev_phase + 2*M_PI/sps;
        if (wrap_to_pi(phase) < -M_PI * margin && wrap_to_pi(prev_phase) > M_PI * margin) {
            symbs[bit_len] = sign(filtered_samps[i]);
            bit_len++;
        }
        if (zc[i-1]){
            error = wrap_to_pi(phase);
            integrator = integrator + error * ki_PLL;
            sps = guess_sps + error*kp_PLL + integrator;
        }
        prev_phase = phase;
    }

    return bit_len;
}