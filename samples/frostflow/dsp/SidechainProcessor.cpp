#include "SidechainProcessor.h"
#include <cmath>
#include <algorithm> // for std::clamp

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace samples {
namespace frostflow {

void SidechainProcessor::init() {
    envelope = 0.0f;
    lfoPhase = 0.0f;
}

float SidechainProcessor::process(float input, float depth, Mode mode, float paramTarget, bool& triggerOutput) {
    float modulation = 0.0f;
    triggerOutput = false;

    // 1. Signal conditioning (rectify for envelope?)
    float absInput = std::abs(input);

    switch (mode) {
        case ENV:
            // Simple AR Envelope Follower
            if (absInput > envelope) {
                envelope += attack * (absInput - envelope);
            } else {
                envelope += release * (absInput - envelope);
            }
            modulation = envelope;
            break;

        case GATE:
            // Schmidt Trigger or simple threshold
            if (!gateState && input > gateThreshold) {
                gateState = true;
                triggerOutput = true;
                modulation = 1.0f; // High
            } else if (gateState && input < (gateThreshold - 0.1f)) {
                gateState = false;
                modulation = 0.0f; // Low
            } else {
                modulation = gateState ? 1.0f : 0.0f;
            }
            break;

        case LFO:
            // Rate modulated by input? Or Input is LFO?
            // Spec says: "LFO: Low-rate CV modulates Dry/Wet continuously"
            // This might mean Internal LFO modulated by Input, or Input IS the LFO source?
            // "SC IN jack (audio/CV input)"
            // "SC MODE switch... LFO: Low-rate CV modulates Dry/Wet continuously"
            // Usually this means the Sidechain Input IS the LFO.
            // But if nothing is plugged in? Maybe an internal LFO?
            // Let's assume Input IS the modulation source if plugged in.
            // If internal LFO is needed, we'd generate it.
            // Let's assume Input is the LFO signal.
            modulation = input * 0.2f; // Scale raw CV (usually -5/5V -> -1/1)
            break;
    }

    // Apply Modulation to Parameter
    // Target is 0..1
    // Modulation is 0..1 (Env/Gate) or Bipolar (LFO)
    // Depth is -1..1
    
    float delta = modulation * depth;
    
    float val = paramTarget + delta;
    if (val < 0.0f) val = 0.0f;
    if (val > 1.0f) val = 1.0f;
    return val;
}

}
}
