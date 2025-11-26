#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/dsp/SpectralFreeze.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SpectralFreeze_Sustain) {
    SpectralFreeze freeze;
    freeze.init();
    
    float outL, outR;
    
    // 1. Feed input signal (Sine) for a while to fill buffers
    // 1024 points. 48k. Sine 440.
    for (size_t i=0; i<2000; ++i) {
        float samp = std::sin(2.0f * M_PI * 440.0f * i / 48000.0f);
        freeze.process(samp, samp, outL, outR, false);
    }
    
    // 2. Freeze
    freeze.process(0.0f, 0.0f, outL, outR, true);
    
    // 3. Feed Silence while Frozen
    // Output should be non-zero (Sustain)
    float maxOut = 0.0f;
    for (size_t i=0; i<5000; ++i) {
        freeze.process(0.0f, 0.0f, outL, outR, true);
        if (std::abs(outL) > maxOut) maxOut = std::abs(outL);
    }
    
    // We expect sound to persist.
    // Note: Depending on window overlaps, it might fluctuate, but max should be significant.
    erb_ASSERT_MSG(maxOut > 0.1f, "Frozen signal did not sustain silence input");
}

erb_TEST_CASE(FrostFlow, SpectralFreeze_Clear) {
    SpectralFreeze freeze;
    freeze.init();
    
    float outL, outR;
    
    // 1. Feed input, Freeze
    for (size_t i=0; i<2000; ++i) {
        float samp = 0.5f;
        freeze.process(samp, samp, outL, outR, false);
    }
    freeze.process(0.0f, 0.0f, outL, outR, true);
    
    // 2. Unfreeze and feed silence
    // Output should decay to zero (latency of window/overlap)
    // FFT Size 1024. OLA latency ~ 1024 samples.
    
    for (size_t i=0; i<4096; ++i) {
        freeze.process(0.0f, 0.0f, outL, outR, false);
    }
    
    // Check for silence
    float maxOut = 0.0f;
    for (size_t i=0; i<1024; ++i) {
        freeze.process(0.0f, 0.0f, outL, outR, false);
        if (std::abs(outL) > maxOut) maxOut = std::abs(outL);
    }
    
    erb_ASSERT_FLOAT_EQ(maxOut, 0.0f, 0.001f, "Signal did not clear after unfreezing");
}

erb_TEST_CASE(FrostFlow, SpectralFreeze_StereoIndependence) {
    SpectralFreeze freeze;
    freeze.init();
    
    float outL, outR;
    
    // Feed Left only
    for (size_t i=0; i<2000; ++i) {
        freeze.process(0.5f, 0.0f, outL, outR, false);
    }
    
    // Check output
    // L should be active, R silent
    float maxL = 0.0f;
    float maxR = 0.0f;
    
    // Wait for latency
    for (size_t i=0; i<100; ++i) {
        freeze.process(0.5f, 0.0f, outL, outR, false);
        if (std::abs(outL) > maxL) maxL = std::abs(outL);
        if (std::abs(outR) > maxR) maxR = std::abs(outR);
    }
    
    erb_ASSERT_MSG(maxL > 0.1f, "Left channel silent");
    erb_ASSERT_FLOAT_EQ(maxR, 0.0f, 0.001f, "Right channel bleed");
}
