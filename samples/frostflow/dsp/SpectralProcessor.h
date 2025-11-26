#pragma once

#include <vector>
#include "RngWrapper.h"

namespace samples {
namespace frostflow {

class SpectralProcessor {
public:
    void init();
    
    /**
     * Apply spectral processing to magnitude and phase arrays.
     * @param mag Magnitude array (size FFT_SIZE)
     * @param phase Phase array (size FFT_SIZE)
     * @param fftSize Size of the arrays (must be 1024)
     */
    void process(float* mag, float* phase, size_t fftSize);

    // Parameters
    void setTilt(float t) { tilt = t; }   // -1.0 (Low) to 1.0 (High)
    void setDecay(float d) { decay = d; } // 0.0 (Infinite) to 1.0 (Quick)
    
    void setMotionAmount(float m) { motionAmount = m; }

private:
    float tilt = 0.0f;
    float decay = 0.0f;
    float motionAmount = 0.2f; // Default small amount

    RngWrapper rng;
    
    // Phase jitter state per bin for motion
    std::vector<float> phaseOffsets;
};

}
}
