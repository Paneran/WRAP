#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "tools.h"
#include "constants.h"
#include "receiver.h"

void normalize(const uint16_t * samples, float * norm_samples);
void costas_loop(float * samples, float * samples_d, params_r * params);
int timing_recovery(float * filtered_samps, int * symbs, params_r * params);

/*
Takes samples and turns them into symbols.
samples: samples from ADC
symbs: symbs at the output
params: parameters that need to be stored intermediately
returns length of symbs array. symbs array must be 
allocated for longer than samples/sps + some margin
*/

int demodulate(const uint16_t * samples, int * symbs, params_r * params) {
    float norm_samples[NUM_SAMPLES];
    normalize(samples, norm_samples);
    
    if (DEBUG) {
        FILE *fpt1 = fopen("data/norm_wave.csv", "w+");
        for (int i = 0; i < NUM_SAMPLES; i++) {
            fprintf(fpt1, "%f,", norm_samples[i]);
        }
        fclose(fpt1);
    }

    // Costas Loop
    float samples_d[NUM_SAMPLES] = {0, 0, 0, 0, 0, 0};
    costas_loop(norm_samples, samples_d, params);

    if (DEBUG) {
        FILE *fpt3 = fopen("data/costa_samples.csv", "w+");
        for (int i = 0; i < NUM_SAMPLES; i++) {
            fprintf(fpt3, "%f,", samples_d[i]);
        }
        fclose(fpt3);
    }
    
    // filter w SRRC
    float filtered_samps[1<<((int)ceil(log2(NUM_SAMPLES)))];
    conv(samples_d, RRC, filtered_samps, NUM_SAMPLES, RRC_LEN);
    // readjust window
    float shift = RRC_LEN/2. - 0.5;
    int k;
    for (int i = shift ; i < NUM_SAMPLES+RRC_LEN-1-shift; i++) {
        k = i - shift;
        filtered_samps[k] = filtered_samps[i];
    }

    if(DEBUG) {
        FILE *fpt2 = fopen("data/SRRC.csv", "w+");
        for (int i = 0; i < NUM_SAMPLES; i++) {
            fprintf(fpt2, "%f,", filtered_samps[i]);
        }
        fclose(fpt2);
    }


    // timing recovery
    int bit_len = timing_recovery(filtered_samps, symbs, params);

    if(DEBUG) {
        FILE *fpt4 = fopen("data/sampled.csv", "w+");
        for (int i = 0; i < N; i++) {
            fprintf(fpt4, "%i,", symbs[i]);
        }
        fclose(fpt4);
    }
    return bit_len;
}


uint8_t find_packet(int * symbs, uint8_t * bits, const int num_symbs) {
    // take cross correlation
    int xcorr_out[CORRELATION_BUFFER];
    uint8_t packet_found = 0;
    xcorr(key, symbs, xcorr_out, 15, num_symbs);

    // find packet
    int shift = 0;
    for (int i = 0; i < num_symbs-30; i++) {
        if (abs(xcorr_out[i] + xcorr_out[i+15] + xcorr_out[i+30] ) > 39) {
            shift = i+31;
            if (xcorr_out[i] < 0) {
                for (int i = 0; i < N-45; i++) {
                    symbs[shift + i] = symbs[shift+ i]*-1;
                }
            }
            packet_found = 1;
            break;
        }
    }
    
    if (!packet_found) 
        return 0;

    // convert symbols to bits
    for (int i = 0; i < N-45; i++) {
        bits[i] = (symbs[shift+i]+1)*0.5;
    }
    return 1;
}

void normalize(const uint16_t * samples, float * norm_samples) {
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
    // divide by 60 arbitrary, just done to get to an ampltiude I used to tune gain values
    for (int i = 0; i < NUM_SAMPLES; i++) {
        norm_samples[i] = (samples[i]-mean)/var*3;
    }
}

void costas_loop(float * norm_samples, float * samples_d, params_r * params) {
    float phase = params->CL_phase;
    float inph[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float quad[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float inph_[ORDER+1] = {0, 0, 0, 0, 0, 0};
    float quad_[ORDER+1] = {0, 0, 0, 0, 0, 0};
    double error = 0;
    double integrator = 0;

    float kp = 8.5;
    float ki = 0.1;
    float dt = FC/FS;

    for (int i = ORDER; i < NUM_SAMPLES+ORDER; i++) {
        // define t from microcontroller
        int k = i - ORDER;
        inph_[ORDER] = norm_samples[k]*2*cos(2*M_PI*dt*k + phase);
        quad_[ORDER] = norm_samples[k]*-2*sin(2*M_PI*dt*k + phase);

        filter(inph_,lp, inph, ORDER+1,ORDER+1);
        filter(quad_,lp, quad, ORDER+1,ORDER+1);

        samples_d[k] = inph[ORDER];

        error = inph[ORDER] * quad[ORDER];
        integrator = integrator + ki*error;
        phase = phase + kp*error + integrator;
        
        // shift the values of inph_ and quad_
        for (int j = 1; j < ORDER+1; j++) {
            inph_[j-1] = inph_[j];
            quad_[j-1] = quad_[j];
        }
    }
}

int timing_recovery(float * filtered_samps, int * symbs, params_r * params) {
    const float kp_PLL = 0.01;
    const float ki_PLL = 0.05;
    const float margin = 0.75;

    float sps = params->sps;
    int integrator = 0, error = 0, bit_len = 0;
    float prev_phase = params->TR_phase, phase;
    
    float zc[NUM_SAMPLES-1];

    // calculate zero crossings
    for (int i = 0; i < NUM_SAMPLES-1; i++) {
        int temp = (filtered_samps[i+1] * filtered_samps[i])/abs((filtered_samps[i+1] * filtered_samps[i]));
        zc[i] = !(temp+1);
    }

    // timing recovery
    for (int i = 1; i < NUM_SAMPLES; i++) {
        phase = prev_phase + 2*M_PI/sps;
        phase = wrap_to_pi(phase);
        if (phase < -M_PI * margin && prev_phase > M_PI * margin) {
            symbs[bit_len] = (int)(filtered_samps[i]/fabs(filtered_samps[i]));
            bit_len++;
        }
        if (bit_len==N) {
            break;
        }
        if (zc[i-1]){
            error = phase;
            integrator = integrator + error * ki_PLL;
            sps = SPS + error*kp_PLL + integrator;
        }
        prev_phase = phase;
    }
    params->sps = sps;
    params->TR_phase = phase;
    return bit_len;
}