#include "BlendEngine.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace samples {
namespace frostflow {

void BlendEngine::init() {
}

void BlendEngine::equalPowerFade(float x, float& outA, float& outB) {
    // Standard Equal Power Crossfade
    // Using sin/cos law: gain^2 + gain^2 = 1
    float angle = x * (M_PI * 0.5f);
    outA = std::cos(angle);
    outB = std::sin(angle);
}

void BlendEngine::process(float dryWet, float abMorph, float& gainDry, float& gainA, float& gainB) {
    // 1. Calculate Mix between Dry and Wet(Layers)
    float mixDry = 0.0f;
    float mixWet = 0.0f;
    equalPowerFade(dryWet, mixDry, mixWet);
    
    // 2. Calculate Mix between Layer A and Layer B
    float mixA = 0.0f;
    float mixB = 0.0f;
    equalPowerFade(abMorph, mixA, mixB);

    // 3. Combine
    gainDry = mixDry;
    gainA = mixWet * mixA;
    gainB = mixWet * mixB;
}

}
}
