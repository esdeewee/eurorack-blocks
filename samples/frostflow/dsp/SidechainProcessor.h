#pragma once

#include <cmath>

namespace samples {
namespace frostflow {

class SidechainProcessor {
public:
    enum Mode {
        ENV,
        GATE,
        LFO
    };

    void init();

    /**
     * Process sidechain input and modulation.
     * @param input Audio/CV input
     * @param depth -1.0 to 1.0 modulation depth
     * @param mode ENV, GATE, or LFO
     * @param paramTarget The parameter to modulate (e.g. Dry/Wet)
     * @param triggerOutput Output for gate triggers (true if triggered)
     * @return Modulated parameter value
     */
    float process(float input, float depth, Mode mode, float paramTarget, bool& triggerOutput);

private:
    // Envelope Follower
    float envelope = 0.0f;
    float attack = 0.01f; // coeff
    float release = 0.001f; // coeff

    // Gate Detector
    bool gateState = false;
    float gateThreshold = 0.5f;

    // LFO (Simple Sine for now)
    float lfoPhase = 0.0f;
    float lfoRate = 0.0001f; // per sample
};

}
}
