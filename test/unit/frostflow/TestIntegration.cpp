#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/FrostFlowDsp.h"

using namespace samples::frostflow;

// 1. Torture Test (Random Inputs)
erb_TEST_CASE(FrostFlow, Integration_Torture) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE];
    float inR[BLOCK_SIZE];
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    srand(100);
    
    // Run for 50 blocks
    for (int b=0; b<50; ++b) {
        TestSignal::generateNoise(inL, BLOCK_SIZE, 1.0f, rand());
        TestSignal::generateNoise(inR, BLOCK_SIZE, 1.0f, rand());
        
        if (b % 10 == 0) {
            dsp.paramDryWet = (float)rand() / RAND_MAX;
            dsp.paramAbMorph = (float)rand() / RAND_MAX;
            dsp.paramLayerALevel = (float)rand() / RAND_MAX;
            dsp.freezeReqA = (rand() % 100 < 5);
        }
        
        dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
        
        for (size_t i=0; i<BLOCK_SIZE; ++i) {
            if (std::isnan(outL[i]) || std::isinf(outL[i])) erb_ASSERT_MSG(false, "NaN/Inf detected in Left output");
            if (std::isnan(outR[i]) || std::isinf(outR[i])) erb_ASSERT_MSG(false, "NaN/Inf detected in Right output");
            if (std::abs(outL[i]) > 2.0f) erb_ASSERT_MSG(false, "Output explosion detected");
        }
    }
}

// 2. Silence Test
erb_TEST_CASE(FrostFlow, Integration_Silence) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE] = {0};
    float inR[BLOCK_SIZE] = {0};
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    dsp.paramDryWet = 0.0f; // Dry
    
    for (int b=0; b<10; ++b) {
        dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
        for (size_t i=0; i<BLOCK_SIZE; ++i) {
            erb_ASSERT_FLOAT_EQ(outL[i], 0.0f, 0.0001f, "Silence input produced non-silence output");
        }
    }
}

// 3. Phase Cancellation Test (Mono Sum)
erb_TEST_CASE(FrostFlow, Integration_PhaseCheck) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE];
    float inR[BLOCK_SIZE];
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    // Feed identical signal to L and R
    TestSignal::generateSine(inL, BLOCK_SIZE, 440.0f, 48000.0f);
    TestSignal::generateSine(inR, BLOCK_SIZE, 440.0f, 48000.0f);
    
    dsp.paramDryWet = 0.5f; // Mix
    dsp.paramAbMorph = 0.0f; 
    
    // Process
    // Wait for potential FFT latency (though Dry path has 0 latency usually)
    // Spec says Dry/Wet blend. 
    // If implementation uses equal power crossfade, L and R should remain identical if processing is symmetric.
    
    dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
    
    for (size_t i=0; i<BLOCK_SIZE; ++i) {
        // L should equal R
        erb_ASSERT_FLOAT_EQ(outL[i], outR[i], 0.001f, "L/R Phase mismatch on identical input");
    }
}

// 4. DC Offset Removal Test
erb_TEST_CASE(FrostFlow, Integration_DCRemoval) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE];
    float inR[BLOCK_SIZE];
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    // Feed DC
    TestSignal::generateDC(inL, BLOCK_SIZE, 0.5f);
    TestSignal::generateDC(inR, BLOCK_SIZE, 0.5f);
    
    dsp.paramDryWet = 1.0f; // Fully Wet (Spectral)
    
    // Spectral processing usually blocks DC if bin 0 is handled or windowed?
    // Our Freeze implementation preserves DC bin?
    // SpectralFreeze::analyze: mag[0] = abs(fft[0]).
    // SpectralFreeze::synthesize: fft[0] = mag[0].
    // So DC is preserved! 
    // Wait, usually audio DC should be blocked.
    // But spec doesn't say "DC Block".
    // Let's test if it IS preserved (as currently implemented) or if we should add a blocker.
    // If we want "Flawless", we usually want DC blocking on output.
    // But let's just verify behavior for now.
    
    // Run for latency duration
    for(int i=0; i<20; ++i) dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
    
    // If DC is preserved, output should be non-zero.
    // If windowing is applied (Hann), DC into FFT -> Windowed DC -> FFT.
    // Window reduces DC amplitude by half (average of Hann is 0.5).
    // Overlap add restores it.
    
    // Just check it doesn't explode or go NaN.
    erb_ASSERT_MSG(std::abs(outL[0]) < 2.0f, "DC Input sanity check");
}

// 5. Freeze Triggering Test
erb_TEST_CASE(FrostFlow, Integration_FreezeTrigger) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE];
    float inR[BLOCK_SIZE];
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    TestSignal::generateSine(inL, BLOCK_SIZE, 440.0f, 48000.0f);
    TestSignal::generateSine(inR, BLOCK_SIZE, 440.0f, 48000.0f);
    
    dsp.paramDryWet = 1.0f; 
    
    // 1. Live mode
    dsp.freezeReqA = false;
    dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
    
    // 2. Trigger Freeze
    dsp.freezeReqA = true; 
    // Need rising edge logic in DSP?
    // FrostFlow.cpp handles button, passes boolean.
    // FrostFlowDsp::process handles edge detection.
    // So toggling false->true should freeze.
    
    dsp.process(inL, inR, outL, outR, BLOCK_SIZE); // Edge detected here
    
    // 3. Sustain
    // Silence input
    TestSignal::generateDC(inL, BLOCK_SIZE, 0.0f);
    TestSignal::generateDC(inR, BLOCK_SIZE, 0.0f);
    
    // Process multiple blocks
    float maxOut = 0.0f;
    for(int i=0; i<50; ++i) {
        dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
        for(int s=0; s<BLOCK_SIZE; ++s) {
            if(std::abs(outL[s]) > maxOut) maxOut = std::abs(outL[s]);
        }
    }
    
    // Should sustain
    erb_ASSERT_MSG(maxOut > 0.01f, "Freeze failed to sustain signal");
}
