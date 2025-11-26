#pragma once
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace samples {
namespace frostflow {

class Limiter {
public:
    void init() {
        currentGain = 1.0f;
    }

    // Simple lookahead limiter or soft clipper
    // For Eurorack context, a soft clipper/saturator is often more musical than a hard brickwall.
    // Let's implement a soft knee limiter.
    void process(float* audio, size_t size) {
        const float threshold = 0.95f;
        const float release = 0.9995f; // Slow release
        
        for (size_t i = 0; i < size; ++i) {
            float absIn = std::abs(audio[i]);
            
            if (absIn > threshold) {
                // Instant attack
                float targetGain = threshold / absIn;
                if (targetGain < currentGain) {
                    currentGain = targetGain;
                }
            } else {
                // Release
                currentGain = currentGain * release + 1.0f * (1.0f - release);
            }
            
            // Apply
            audio[i] *= currentGain;
            
            // Safety Hard Clip just in case
            if (audio[i] > 1.0f) audio[i] = 1.0f;
            if (audio[i] < -1.0f) audio[i] = -1.0f;
        }
    }

private:
    float currentGain = 1.0f;
};

}
}
