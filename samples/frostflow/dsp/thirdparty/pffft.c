/*
   PFFFT : a Pretty Fast FFT.
   (Simulated version for Unit Testing if real source missing)
*/

#include "../../submodules/vcv-rack-sdk/dep/include/pffft.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct PFFFT_Setup {
    int N;
    pffft_transform_t transform;
};

PFFFT_Setup *pffft_new_setup(int N, pffft_transform_t transform) {
    PFFFT_Setup *s = (PFFFT_Setup*)malloc(sizeof(PFFFT_Setup));
    s->N = N;
    s->transform = transform;
    return s;
}

void pffft_destroy_setup(PFFFT_Setup *s) {
    free(s);
}

void *pffft_aligned_malloc(size_t nb_bytes) {
    return malloc(nb_bytes);
}

void pffft_aligned_free(void *p) {
    free(p);
}

int pffft_simd_size() { return 1; }

void pffft_transform(PFFFT_Setup *setup, const float *input, float *output, float *work, pffft_direction_t direction) {
    // Not implemented
}

void pffft_transform_ordered(PFFFT_Setup *setup, const float *input, float *output, float *work, pffft_direction_t direction) {
    int N = setup->N;
    if (setup->transform != PFFFT_REAL) return; 

    if (direction == PFFFT_FORWARD) {
        // DC
        double sum_re = 0;
        for (int n = 0; n < N; ++n) sum_re += input[n];
        output[0] = (float)sum_re;
        
        // Nyquist (k=N/2)
        sum_re = 0;
        for (int n = 0; n < N; ++n) {
            sum_re += input[n] * ((n % 2) ? -1.0 : 1.0);
        }
        output[1] = (float)sum_re;

        for (int k = 1; k < N / 2; ++k) {
            sum_re = 0;
            double sum_im = 0;
            for (int n = 0; n < N; ++n) {
                double angle = -2.0 * M_PI * k * n / N;
                double c = cos(angle);
                double s = sin(angle);
                sum_re += input[n] * c;
                sum_im += input[n] * s;
            }
            output[2 * k] = (float)sum_re;
            output[2 * k + 1] = (float)sum_im;
        }
    } else {
        // Inverse
        for (int n = 0; n < N; ++n) {
            double sum = 0;
            
            // DC
            sum += input[0];
            
            // Nyquist
            sum += input[1] * ((n % 2) ? -1.0 : 1.0);
            
            for (int k = 1; k < N / 2; ++k) {
                float re = input[2 * k];
                float im = input[2 * k + 1];
                
                double angle = 2.0 * M_PI * k * n / N;
                // Add X[k] + X[N-k] = 2 * Real(X[k] * exp(jwn))
                sum += 2.0 * (re * cos(angle) - im * sin(angle));
            }
            output[n] = (float)sum; 
        }
    }
}

void pffft_zreorder(PFFFT_Setup *setup, const float *input, float *output, pffft_direction_t direction) {
}

void pffft_zconvolve_accumulate(PFFFT_Setup *setup, const float *dft_a, const float *dft_b, float *dft_ab, float scaling) {
}
