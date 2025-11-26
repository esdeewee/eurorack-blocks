#include "FrostFlowDsp.h"

namespace samples {
namespace frostflow {

void FrostFlowDsp::init() {
    layerA.init();
    layerB.init();
    blendEngine.init();
    sidechainProcessor.init();
    clockQuantizerA.init();
    clockQuantizerB.init();
    
    isFrozenA = false;
    isFrozenB = false;
}

void FrostFlowDsp::process(const float* in_l, const float* in_r, float* out_l, float* out_r, size_t size) {
    
    // Sidechain Processing
    bool triggerOut = false;
    float modulatedDryWet = sidechainProcessor.process(sidechainInput, paramSidechainDepth, sidechainMode, paramDryWet, triggerOut);
    
    // Edge Detection for Freeze Requests
    bool reqA = freezeReqA && !lastFreezeReqA;
    lastFreezeReqA = freezeReqA;
    
    bool reqB = freezeReqB && !lastFreezeReqB;
    lastFreezeReqB = freezeReqB;
    
    // Clock Quantizer Processing
    bool triggerA = clockQuantizerA.process(clockGate, clockDiv, sampleRate, reqA);
    if (triggerA) isFrozenA = !isFrozenA;
    
    bool triggerB = clockQuantizerB.process(clockGate, clockDiv, sampleRate, reqB);
    if (triggerB) isFrozenB = !isFrozenB;

    // Blend Processing
    float gainDry, gainA, gainB;
    blendEngine.process(modulatedDryWet, paramAbMorph, gainDry, gainA, gainB);
    
    // Apply Layer Levels
    gainA *= paramLayerALevel;
    gainB *= paramLayerBLevel;

    // Audio Processing Loop
    for (size_t i = 0; i < size; ++i) {
        float inL = in_l[i];
        float inR = in_r[i];
        
        float outAL = 0.f, outAR = 0.f;
        layerA.process(inL, inR, outAL, outAR, isFrozenA);

        float outBL = 0.f, outBR = 0.f;
        layerB.process(inL, inR, outBL, outBR, isFrozenB);

        // Mix
        out_l[i] = (inL * gainDry) + (outAL * gainA) + (outBL * gainB);
        out_r[i] = (inR * gainDry) + (outAR * gainA) + (outBR * gainB);
    }
}

} // namespace frostflow
} // namespace samples
