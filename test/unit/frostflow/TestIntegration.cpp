#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/FrostFlowDsp.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, Integration_Torture) {
    FrostFlowDsp dsp;
    dsp.init();
    
    const size_t BLOCK_SIZE = 64;
    float inL[BLOCK_SIZE];
    float inR[BLOCK_SIZE];
    float outL[BLOCK_SIZE];
    float outR[BLOCK_SIZE];
    
    srand(100);
    
    // Run for 50 blocks (approx 0.1 seconds)
    // Reduced from 1000 because Mock FFT is O(N^2) and slow
    for (int b=0; b<50; ++b) {
        // 1. Random Audio
        TestSignal::generateNoise(inL, BLOCK_SIZE, 1.0f, rand());
        TestSignal::generateNoise(inR, BLOCK_SIZE, 1.0f, rand());
        
        // 2. Random Params
        if (b % 10 == 0) {
            dsp.paramDryWet = (float)rand() / RAND_MAX;
            dsp.paramAbMorph = (float)rand() / RAND_MAX;
            dsp.paramLayerALevel = (float)rand() / RAND_MAX;
            dsp.freezeReqA = (rand() % 100 < 5); // Occasional freeze
        }
        
        // 3. Process
        dsp.process(inL, inR, outL, outR, BLOCK_SIZE);
        
        // 4. Check Output Sanity (No NaNs, No Inf)
        for (size_t i=0; i<BLOCK_SIZE; ++i) {
            if (std::isnan(outL[i]) || std::isinf(outL[i])) {
                erb_ASSERT_MSG(false, "NaN/Inf detected in Left output");
            }
            if (std::isnan(outR[i]) || std::isinf(outR[i])) {
                erb_ASSERT_MSG(false, "NaN/Inf detected in Right output");
            }
            // Check bounds (Limiter should keep it <= 1.0, but maybe slight overshoot allowed)
            // Let's say < 2.0 is safe 'sanity'.
            if (std::abs(outL[i]) > 2.0f) erb_ASSERT_MSG(false, "Output explosion detected");
        }
    }
}

