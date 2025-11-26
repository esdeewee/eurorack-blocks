#include "SpectralProcessor.h"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace samples {
namespace frostflow {

void SpectralProcessor::init() {
    // Initialize phase offsets
    // We assume 1024 size for now, but better to resize dynamically if needed.
    // For embedded, fixed size is better.
    phaseOffsets.resize(1024, 0.0f);
    rng.setSeed(1234); // Deterministic seed
}

void SpectralProcessor::process(float* mag, float* phase, size_t fftSize) {
    // 1. Decay
    // Decay is 0.0 (Infinite hold) to 1.0 (Quick fade).
    // We multiply magnitudes by (1.0 - decay * factor).
    // If decay is 0, multiplier is 1.0 (no change).
    // If decay is 1, multiplier is small (fast fade).
    // Let's map decay parameter to a per-frame attenuation coefficient.
    // Frame rate ~ 48000 / 256 ~ 187 Hz.
    // Decay coeff: 0.999 (long) to 0.8 (short).
    // Map 0..1 to that range.
    // Actually spec says: "Decay envelope (infinite hold -> quick fade)"
    // If 'frozen', we are holding.
    // So this process() is called continuously on the frozen buffer?
    // If so, we modify the buffer in place.
    // Wait, if we modify in place, the buffer eventually silences.
    // So 'infinite hold' means we don't modify.
    
    float decayCoeff = 1.0f;
    if (decay > 0.01f) {
        // Map 0..1 linear to exponential decay
        // min coeff 0.90 (fast fade), max 1.0
        decayCoeff = 1.0f - (decay * decay * 0.1f); 
    }

    // 2. Tilt
    // Tilt -1 (Low boost, High cut) to +1 (Low cut, High boost)
    // Simple linear slope in dB or magnitude?
    // Let's use a spectral slope.
    // Bin 0 is DC, Bin N/2 is Nyquist.
    // Slope = bin_index / (N/2).
    
    float tiltSlope = tilt; // -1..1

    // 3. Motion (Phase Blur/Jitter)
    // We add random phase offsets to blur the sound.
    // Motion amount controls the speed of phase drift.
    
    for (size_t k = 0; k < fftSize; ++k) {
        // Apply Decay
        mag[k] *= decayCoeff;

        // Apply Tilt
        // Normalized freq 0..1
        float freqNorm = (float)k / (float)(fftSize / 2);
        // Limit to 1.0 max (mirroring for > N/2 logic if we processed full array?
        // PFFFT Ordered: 0=DC, 1=Nyquist.
        // 2k, 2k+1 are Re/Im of bin k.
        // We are passed vectors of Mag/Phase which are size N.
        // Our SpectralFreeze analysis loop filled:
        // mag[0], mag[1]... mag[k]...
        // Indices 0 and 1 are special.
        // k=1..(N/2-1) correspond to mag[2*k].
        
        // Wait, SpectralFreeze::analyze fills the vector differently?
        // Let's double check SpectralFreeze::analyze storage format.
        // "mag[2 * k] = std::sqrt(re * re + im * im);"
        // "mag[2 * k + 1] = mag[2 * k];"
        // So indices correspond to 2*bin_index.
        
        float binIndex = 0.0f;
        if (k == 0) binIndex = 0.0f;
        else if (k == 1) binIndex = (float)(fftSize / 2);
        else binIndex = (float)(k / 2);
        
        float f = binIndex / (float)(fftSize / 2);
        
        // Tilt Logic
        // If tilt > 0 (High boost): multiply by f
        // If tilt < 0 (Low boost): multiply by (1-f)
        // This is aggressive. Let's use a gentle shelf or slope.
        // dB slope: +/- 6dB per octave?
        // Simple linear gain ramp:
        // Gain = 1.0 + tilt * (2*f - 1.0)
        // Range 0..2.
        
        float tiltGain = 1.0f;
        if (std::abs(tilt) > 0.01f) {
            // Linear ramp from (1-tilt) to (1+tilt)
            // tilt +1: 0 to 2
            // tilt -1: 2 to 0
            tiltGain = 1.0f + tilt * (2.0f * f - 1.0f);
            if (tiltGain < 0.0f) tiltGain = 0.0f;
            if (tiltGain > 2.0f) tiltGain = 2.0f;
        }
        mag[k] *= tiltGain;

        // Apply Motion (Phase Jitter)
        // Only apply to phase bins, not magnitudes (though passed in same loop)
        // We have phase[k].
        // DC (0) and Nyquist (1) phases are usually 0 or Pi, fixed.
        // Bins k > 1 are phases.
        
        if (k > 1) {
            // Random walk phase
            // We need persistent state for "Drift".
            // phaseOffsets[k] += rng.nextBiFloat() * motionAmount * 0.1f;
            // phase[k] += phaseOffsets[k];
            
            // For "Instant Musicality", phase blur is good.
            // Just adding random jitter every frame makes it noise?
            // No, slow drift.
            
            float jitter = rng.nextBiFloat() * motionAmount * 0.2f;
            phaseOffsets[k] += jitter;
            
            // Wrap phase
            if (phaseOffsets[k] > M_PI) phaseOffsets[k] -= 2.0f * M_PI;
            if (phaseOffsets[k] < -M_PI) phaseOffsets[k] += 2.0f * M_PI;
            
            phase[k] += phaseOffsets[k];
        }
    }
}

}
}
