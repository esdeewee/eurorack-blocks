#pragma once

#include <cstddef>
#include <vector>

#if defined(erb_TARGET_DAISY)
// #include "arm_math.h" // Phase 2
#else
#include <pffft.h>
#endif

namespace samples {
namespace frostflow {

class SpectralFft {
public:
    static const size_t FFT_SIZE = 1024;

    SpectralFft();
    ~SpectralFft();

    void init();
    
    /**
     * Perform Forward Real FFT
     * @param input Time domain signal (must be FFT_SIZE length)
     * @param output Frequency domain signal (must be FFT_SIZE length)
     * 
     * Output format (Ordered):
     * bin 0: Real(0) (DC)
     * bin 1: Real(N/2) (Nyquist)
     * bin 2k: Real(k)
     * bin 2k+1: Imag(k)
     * for k = 1 .. N/2 - 1
     */
    void forward(const float* input, float* output);

    /**
     * Perform Inverse Real FFT
     * @param input Frequency domain signal (must be FFT_SIZE length)
     * @param output Time domain signal (must be FFT_SIZE length)
     */
    void inverse(const float* input, float* output);

private:
#if !defined(erb_TARGET_DAISY)
    PFFFT_Setup* setup = nullptr;
    float* workBuffer = nullptr;
#endif
};

}
}
