#pragma once

#include <array>
#include "dsp/SpectralFft.h"
#include "dsp/SpectralFreeze.h"
#include "dsp/BlendEngine.h"
#include "dsp/SidechainProcessor.h"
#include "dsp/ClockQuantizer.h"
#include "dsp/SpectralProcessor.h"
#include "dsp/RngWrapper.h"
#include "dsp/Limiter.h"

namespace samples {
namespace frostflow {

class FrostFlowDsp {
public:
    void init();
    // Inputs are const
    void process(const float* in_l, const float* in_r, float* out_l, float* out_r, size_t size);

    // Parameters
    bool freezeReqA = false;
    bool freezeReqB = false;
    
    // State (ReadOnly for UI)
    bool isFrozenA = false;
    bool isFrozenB = false;
    
    float paramDryWet = 0.0f; // 0..1
    float paramAbMorph = 0.0f; // 0..1
    float paramSidechainDepth = 0.0f; // -1..1
    
    // Sidechain Inputs
    float sidechainInput = 0.0f;
    SidechainProcessor::Mode sidechainMode = SidechainProcessor::ENV;

    // Clock Inputs
    bool clockGate = false;
    ClockQuantizer::Division clockDiv = ClockQuantizer::DIV_4;

    // Layer Params
    float paramLayerALevel = 1.0f;
    float paramLayerBLevel = 1.0f;
    
    // Getters for SpectralFreeze layers to set params
    SpectralFreeze& getLayerA() { return layerA; }
    SpectralFreeze& getLayerB() { return layerB; }

private:
    SpectralFreeze layerA;
    SpectralFreeze layerB;
    BlendEngine blendEngine;
    SidechainProcessor sidechainProcessor;
    ClockQuantizer clockQuantizerA;
    ClockQuantizer clockQuantizerB;
    
    float sampleRate = 48000.0f;
    
    // Edge detection state
    bool lastFreezeReqA = false;
    bool lastFreezeReqB = false;
};

} // namespace frostflow
} // namespace samples
