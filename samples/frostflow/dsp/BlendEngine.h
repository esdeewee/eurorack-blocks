#pragma once

#include <cmath>

namespace samples {
namespace frostflow {

class BlendEngine {
public:
    void init();

    /**
     * Calculate blend gains based on Dry/Wet and A/B knobs.
     * @param dryWet 0.0 (Dry) to 1.0 (Wet)
     * @param abMorph 0.0 (A) to 1.0 (B)
     * @param gainDry Output gain for Dry signal
     * @param gainA Output gain for Layer A
     * @param gainB Output gain for Layer B
     */
    void process(float dryWet, float abMorph, float& gainDry, float& gainA, float& gainB);

private:
    // Helper for equal power crossfade
    // x: 0..1
    // outA: cos(x * pi/2)
    // outB: sin(x * pi/2)
    void equalPowerFade(float x, float& outA, float& outB);
};

}
}
