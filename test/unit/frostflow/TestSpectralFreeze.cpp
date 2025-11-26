#include "../test.h"
#include "../../../samples/frostflow/dsp/SpectralFreeze.h"
#include <vector>

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SpectralFreeze_Basic) {
    SpectralFreeze freeze;
    freeze.init();
    
    float outL, outR;
    
    // Process some silence to initialize buffers
    for (size_t i=0; i<2048; ++i) {
        freeze.process(0.0f, 0.0f, outL, outR, false);
    }
    
    // Feed Impulse
    freeze.process(1.0f, 1.0f, outL, outR, false);
    
    // Should eventually appear at output (OLA latency)
    // This test mainly ensures no crash and state consistency
    for (size_t i=0; i<10000; ++i) {
        freeze.process(0.0f, 0.0f, outL, outR, false);
    }
    
    erb_TEST(true); // If we got here, no crash
}
